import os
import math
from PyQt5.QtWidgets import QAction, QTableWidgetItem, QMessageBox
from PyQt5.QtGui import QColor, QIcon
from qgis.core import (
    QgsProject, QgsCoordinateTransform, QgsSpatialIndex, 
    QgsGeometry, QgsWkbTypes
)
from .cadastral_auditor_dialog import CadastralAuditorDialog
from .analyzer import ConfidenceStringMatcher

class CadastralAuditor:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.dlg = CadastralAuditorDialog()
        self.action = None
        self.vector_layer = None # Memory layer for error vectors

    def _setup_vector_layer(self):
        """Creates or clears the memory layer for error vectors."""
        layer_name = "Auditor Error Vectors"
        
        # Find existing layer
        existing_layers = QgsProject.instance().mapLayersByName(layer_name)
        if existing_layers:
            self.vector_layer = existing_layers[0]
            self.vector_layer.dataProvider().truncate()
        else:
            # Create new memory layer
            from qgis.core import QgsVectorLayer, QgsField
            from PyQt5.QtCore import QVariant
            self.vector_layer = QgsVectorLayer("LineString?crs=" + QgsProject.instance().crs().authid(), layer_name, "memory")
            self.vector_layer.dataProvider().addAttributes([
                QgsField("fid", QVariant.Int),
                QgsField("status", QVariant.String)
            ])
            self.vector_layer.updateFields()
            QgsProject.instance().addMapLayer(self.vector_layer)
            
            # Simple Styling
            from qgis.core import QgsSymbol, QgsSingleSymbolRenderer
            symbol = QgsSymbol.defaultSymbol(self.vector_layer.geometryType())
            symbol.setColor(QColor(255, 0, 0)) # Red
            symbol.setWidth(0.5)
            self.vector_layer.setRenderer(QgsSingleSymbolRenderer(symbol))

    def initGui(self):
        self.action = QAction("Cadastral Auditor", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("&Cadastral Auditor", self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        self.iface.removePluginMenu("&Cadastral Auditor", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        self.dlg.populate_layers()
        self.dlg.show()
        
        # Connect signals
        try:
            self.dlg.btn_analyze.clicked.disconnect()
        except TypeError:
            pass 
        self.dlg.btn_analyze.clicked.connect(self.run_analysis)
        
        # Connect table selection (REQ-2026-02-06)
        try:
            self.dlg.table_results.cellClicked.disconnect()
        except TypeError:
            pass
        self.dlg.table_results.cellClicked.connect(self.zoom_to_selected_feature)

    def zoom_to_selected_feature(self, row, column):
        """
        Zooms to the feature corresponding to the clicked row.
        """
        item = self.dlg.table_results.item(row, 1) # Column 1: MatchID (FID)
        if not item:
            return
            
        try:
            fid = int(item.text())
        except ValueError:
            return
            
        cur_layer = self.dlg.cb_layer_current.currentData()
        if not cur_layer:
            return
            
        # Select and Zoom
        cur_layer.selectByIds([fid])
        self.iface.mapCanvas().zoomToSelected(cur_layer)
        self.iface.mapCanvas().refresh()

    def run_analysis(self):
        self.dlg.clear_table()
        self._setup_vector_layer()
        
        # 1. Get Inputs
        cad_layer = self.dlg.cb_layer_cadastral.currentData()
        cur_layer = self.dlg.cb_layer_current.currentData()
        
        if not cad_layer or not cur_layer:
            QMessageBox.warning(self.dlg, "Error", "Layer not selected.")
            return

        tol_min = self.dlg.sb_tol_min.value()
        exclusion_limit = self.dlg.sb_exclusion_limit.value()
        is_precision = self.dlg.rb_mode_precision.isChecked()
        
        # 2. CRS Transform Setup
        xform = None
        if cad_layer.crs() != cur_layer.crs():
            xform = QgsCoordinateTransform(cur_layer.crs(), cad_layer.crs(), QgsProject.instance())
            
        # 3. Build Spatial Index (Cadastral)
        index = QgsSpatialIndex(cad_layer.getFeatures())
        
        # 4. Initialize Matcher
        matcher = ConfidenceStringMatcher(sigma=tol_min)
        
        # 5. Iterate & Process
        feature_count = cur_layer.featureCount()
        self.dlg.progress_bar.setMaximum(feature_count)
        
        row_idx = 0
        for i, cur_feat in enumerate(cur_layer.getFeatures()):
            self.dlg.progress_bar.setValue(i+1)
            
            geom_curr = cur_feat.geometry()
            if not geom_curr: continue
            
            # Apply Transform for search
            search_geom = QgsGeometry(geom_curr)
            if xform:
                search_geom.transform(xform)
                
            # Filter 1: Spatial Search
            bbox = search_geom.boundingBox().buffered(exclusion_limit)
            candidate_ids = index.intersects(bbox)
            
            best_match_feat = None
            min_dist = float('inf')
            
            # Filter 2: Select Best "Shape-Matched" Candidate (Hausdorff Distance)
            # REQ-2026-02-06: Use Hausdorff distance to find the most similar shape
            # instead of just the nearest point. This prevents jumping to perpendicular neighbors.
            
            best_candidate = None
            min_hausdorff = float('inf')
            
            # Pre-filter candidates within a reasonable range (e.g. 2x limit) to reduce cost
            potential_candidates = []
            for cid in candidate_ids:
                cad_feat = cad_layer.getFeature(cid)
                cad_geom = cad_feat.geometry()
                
                # Loose distance check first
                if search_geom.distance(cad_geom) <= exclusion_limit * 2.0:
                    potential_candidates.append(cad_feat)
            
            # Find the single best candidate based on Hausdorff Distance
            for cad_feat in potential_candidates:
                cad_geom = cad_feat.geometry()
                h_dist = search_geom.hausdorffDistance(cad_geom)
                
                if h_dist < min_hausdorff:
                    min_hausdorff = h_dist
                    best_candidate = cad_feat
            
            if best_candidate and min_hausdorff < exclusion_limit * 5.0: # Sanity check
                # Execute Analysis with the LOCKED Best Candidate
                # Pass as a list because process_multi_context expects a list
                result = matcher.process_multi_context(
                    cur_feat, 
                    [best_candidate], 
                    transform=xform,
                    tolerance=tol_min
                )
                
                self._add_result_row(row_idx, cur_feat.id(), result)
                
                # Add Visualization (REQ-2026-02-06)
                if "error_vectors" in result:
                    self._add_error_vector(cur_feat.id(), result["status"], result["error_vectors"])
                elif "error_line" in result:
                    self._add_error_vector(cur_feat.id(), result["status"], result["error_line"])
                
                row_idx += 1
        
        # Finalize display
        if self.vector_layer:
            self.vector_layer.triggerRepaint()
            self.iface.mapCanvas().refresh()

    def _add_error_vector(self, fid, status, geom):
        """Adds a feature to the error vector layer."""
        from qgis.core import QgsFeature
        if not self.vector_layer: return
        
        feat = QgsFeature(self.vector_layer.fields())
        feat.setGeometry(geom)
        feat.setAttributes([fid, status])
        self.vector_layer.dataProvider().addFeatures([feat])
        
    def _add_result_row(self, row, fid, res):
        self.dlg.table_results.insertRow(row)
        
        # Columns: No, MatchID, Topology, Score, Decision, ND Cost, Shift X, Shift Y, Reliability, Note
        items = [
            str(row + 1),
            str(fid),
            res.get("topology", "-"),
            f"{res.get('score', 0):.3f}",
            res.get("status", "-"),
            f"{res.get('nd_cost', 0):.3f}",
            f"{res.get('shift_dx', 0):.3f}",
            f"{res.get('shift_dy', 0):.3f}",
            res.get("reliability", "-"),
            res.get("rec_shift", "-")
        ]
        
        for col, text in enumerate(items):
            item = QTableWidgetItem(text)
            
            # Color coding for Decision
            if col == 4: # Decision
                if "부합" in text:
                    item.setBackground(QColor(200, 255, 200)) # Green
                elif "불부합" in text:
                    item.setBackground(QColor(255, 200, 200)) # Red
                else: 
                    item.setBackground(QColor(255, 255, 200)) # Yellow
            
            if col == 8: # Reliability
                if "높음" in text:
                     item.setBackground(QColor(200, 255, 200))
            
            self.dlg.table_results.setItem(row, col, item)
