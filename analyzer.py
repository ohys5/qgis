import math
from qgis.core import QgsGeometry, QgsPointXY, QgsWkbTypes, QgsSpatialIndex, Qgis

class GeometricMatcher:
    """
    Matches Cadastral feature with Current feature using Hausdorff Distance.
    """
    def __init__(self, cadastral_feat, current_feat):
        self.cadastral_geom = cadastral_feat.geometry()
        self.current_geom = current_feat.geometry()

    def get_hausdorff_distance(self):
        # Returns the maximum distance between the two geometries
        return self.current_geom.hausdorffDistance(self.cadastral_geom)

class BestFitSolver:
    """
    Calculates the optimal shift to align centroids and checks residual error.
    """
    def __init__(self, cadastral_feat, current_feat):
        self.cadastral_geom = cadastral_feat.geometry()
        self.current_geom = current_feat.geometry()
    
    def solve(self):
        # 1. Calculate Centroids
        if self.cadastral_geom.isEmpty() or self.current_geom.isEmpty():
            return 0.0, 0.0, 9999.0

        c_cad = self.cadastral_geom.centroid().asPoint()
        c_cur = self.current_geom.centroid().asPoint()
        
        # 2. Calculate Shift Vector (dx, dy)
        dx = c_cad.x() - c_cur.x()
        dy = c_cad.y() - c_cur.y()
        
        # 3. Create shifted geometry for testing
        shifted_geom = QgsGeometry(self.current_geom)
        shifted_geom.translate(dx, dy)
        
        # 4. Calculate Residual Error (Hausdorff distance after shift)
        residual = shifted_geom.hausdorffDistance(self.cadastral_geom)
        
        return dx, dy, residual

def SegmentAuditor(cadastral_feat, current_feat, tolerance=0.1):
    """
    Analyzes each segment of the current feature against the cadastral feature.
    Returns a summary string of the analysis.
    """
    current_geom = current_feat.geometry()
    cadastral_geom = cadastral_feat.geometry()
    
    # Convert to vertices to simulate segment analysis
    vertices = current_geom.vertices()
    
    issues = []
    for v in vertices:
        point_geom = QgsGeometry.fromPointXY(QgsPointXY(v))
        dist = point_geom.distance(cadastral_geom)
        
        if dist > tolerance:
            # Check if it is inside (gap) or outside (protrusion)
            if cadastral_geom.contains(point_geom):
                issues.append("이격")
            else:
                issues.append("돌출")
                
    return ", ".join(set(issues)) if issues else "부합"

class SmartGeometryComparator:
    def __init__(self, cadastral_feat, current_features):
        self.cadastral_feat = cadastral_feat
        self.current_features = current_features
        self.merged_geom = None
        self.mode = "Line Mode"

    def process(self):
        # 1. Combine and Merge Lines
        geoms = []
        for f in self.current_features:
            g = f.geometry()
            if g and not g.isEmpty():
                geoms.append(g)
        
        if not geoms:
            self.merged_geom = QgsGeometry()
            print("DEBUG: [Analyzer] No valid geometries in current_features")
        else:
            # Initialize with the first geometry
            combined = QgsGeometry(geoms[0])
            for g in geoms[1:]:
                combined = combined.combine(g)
            
            self.merged_geom = combined.mergeLines()
            
            # Fallback
            if not self.merged_geom or self.merged_geom.isEmpty():
                print("DEBUG: [Analyzer] mergeLines() returned empty, using combined geometry")
                self.merged_geom = combined
        
        print(f"DEBUG: [Analyzer] Merged Geom: {self.merged_geom.type()} (Empty: {self.merged_geom.isEmpty()})")
        
        # 2. Auto-Close Detection (5cm)
        if self.merged_geom.type() == QgsWkbTypes.LineGeometry:
            vertices = [v for v in self.merged_geom.vertices()]
            if vertices:
                start = vertices[0]
                end = vertices[-1]
                if start.distance(end) <= 0.05:
                    # Convert to Polygon
                    points = [QgsPointXY(v.x(), v.y()) for v in vertices]
                    if start.distance(end) > 0:
                        points.append(QgsPointXY(start.x(), start.y()))
                    
                    poly_geom = QgsGeometry.fromPolygonXY([points])
                    if poly_geom.isGeosValid():
                        self.merged_geom = poly_geom
                        self.mode = "Polygon Mode"

        result = {"mode": self.mode}
        cad_geom = self.cadastral_feat.geometry()

        if self.mode == "Polygon Mode":
            area_cur = self.merged_geom.area()
            area_cad = cad_geom.area()
            
            intersection = self.merged_geom.intersection(cad_geom)
            overlap_ratio = 0.0
            if area_cad > 0:
                overlap_ratio = (intersection.area() / area_cad) * 100.0
            
            result.update({
                "area_current": area_cur,
                "area_cadastral": area_cad,
                "overlap_ratio": overlap_ratio
            })
        else:
            # Line Mode: Explode and Find Target Segment
            segments = []
            
            print(f"DEBUG: [Analyzer] Cadastral Geom Type: {cad_geom.type()}, WKB: {cad_geom.wkbType()}")
            # Fix: Handle both Polygon and Line geometries to prevent TypeError
            if cad_geom.type() == QgsWkbTypes.PolygonGeometry:
                polys = cad_geom.asMultiPolygon() if cad_geom.isMultipart() else [cad_geom.asPolygon()]
                for poly in polys:
                    for ring in poly:
                        for i in range(len(ring) - 1):
                            segments.append(QgsGeometry.fromPolylineXY([ring[i], ring[i+1]]))
            elif cad_geom.type() == QgsWkbTypes.LineGeometry:
                lines = cad_geom.asMultiPolyline() if cad_geom.isMultipart() else [cad_geom.asPolyline()]
                for line in lines:
                    for i in range(len(line) - 1):
                        segments.append(QgsGeometry.fromPolylineXY([line[i], line[i+1]]))
            
            print(f"DEBUG: [Analyzer] Extracted segments: {len(segments)}")
            
            # Find closest segment (Target Segment)
            target_segment = min(segments, key=lambda s: self.merged_geom.hausdorffDistance(s), default=None)
            
            # Calculate Perpendicular Distance (Centroid to Segment)
            centroid = self.merged_geom.centroid()
            error_line = None
            perp_dist = -1.0
            
            if target_segment and not centroid.isEmpty():
                print(f"DEBUG: [Analyzer] Target Segment: {target_segment.asWkt()}")
                print(f"DEBUG: [Analyzer] Centroid: {centroid.asWkt()}")
                perp_dist = target_segment.distance(centroid)
                closest_pt = target_segment.nearestPoint(centroid)
                if not closest_pt.isEmpty():
                    error_line = QgsGeometry.fromPolylineXY([centroid.asPoint(), closest_pt.asPoint()])
                    print(f"DEBUG: [Analyzer] Error line created: {error_line.asWkt()}")
                else:
                    print("DEBUG: [Analyzer] Closest point is empty")
            else:
                print("DEBUG: [Analyzer] Fallback to Hausdorff (No target segment or empty centroid)")
                perp_dist = self.merged_geom.hausdorffDistance(cad_geom)
            
            print(f"DEBUG: [Analyzer] Calculated perp_dist: {perp_dist}")
            
            result.update({
                "max_discrepancy": perp_dist,
                "error_line": error_line
            })
            
        return result

class ParcelBasedAuditor:
    """
    Analyzes survey lines by intersecting them with cadastral parcels (Polygon).
    Performs Reference Matching and calculates Average Perpendicular Distance.
    """
    def __init__(self, cadastral_features, survey_feature):
        self.cadastral_features = cadastral_features
        self.survey_feature = survey_feature

    def run(self):
        results = []
        survey_geom = self.survey_feature.geometry()
        
        for cad_feat in self.cadastral_features:
            cad_geom = cad_feat.geometry()
            
            if not cad_geom.intersects(survey_geom):
                continue
                
            intersection = survey_geom.intersection(cad_geom)
            if intersection.isEmpty():
                continue
                
            segments = intersection.asMultiPolyline() if intersection.isMultipart() else [intersection.asPolyline()]
            
            for seg_points in segments:
                if len(seg_points) < 2: continue
                seg_geom = QgsGeometry.fromPolylineXY(seg_points)
                if seg_geom.length() == 0: continue
                
                ref_info = self.find_reference_line(seg_geom, cad_geom)
                
                if ref_info:
                    avg_dist = self.calculate_average_distance(seg_geom, ref_info['geometry'])
                    
                    pnu = str(cad_feat.id())
                    for field in ['jibun', 'pnu', 'JIBUN', 'PNU']:
                        if cad_feat.fieldNameIndex(field) != -1:
                            pnu = str(cad_feat[field])
                            break
                    
                    results.append({
                        "pnu": pnu,
                        "direction": ref_info['direction'],
                        "avg_error": avg_dist,
                        "max_error": seg_geom.hausdorffDistance(ref_info['geometry'])
                    })
        return results

    def find_reference_line(self, survey_segment, parcel_geom):
        boundaries = []
        polys = parcel_geom.asMultiPolygon() if parcel_geom.isMultipart() else [parcel_geom.asPolygon()]
        
        for poly in polys:
            for ring in poly:
                for i in range(len(ring) - 1):
                    boundaries.append(QgsGeometry.fromPolylineXY([ring[i], ring[i+1]]))
        
        survey_angle = self._get_angle(survey_segment)
        parcel_center = parcel_geom.centroid().asPoint()
        
        best_ref = None
        min_dist = float('inf')
        best_dir = "Unknown"
        
        for boundary in boundaries:
            b_angle = self._get_angle(boundary)
            diff = abs(survey_angle - b_angle)
            if diff > 180: diff = 360 - diff
            
            if diff < 30 or abs(diff - 180) < 30:
                dist = boundary.distance(survey_segment)
                if dist < min_dist:
                    min_dist = dist
                    best_ref = boundary
                    
                    mid = boundary.interpolate(boundary.length()/2).asPoint()
                    dx = mid.x() - parcel_center.x()
                    dy = mid.y() - parcel_center.y()
                    if abs(dx) > abs(dy):
                        best_dir = "동측" if dx > 0 else "서측"
                    else:
                        best_dir = "북측" if dy > 0 else "남측"
                        
        return {'geometry': best_ref, 'direction': best_dir} if best_ref else None

    def _get_angle(self, line):
        if line.isEmpty(): return 0
        p1 = line.vertexAt(0)
        p2 = line.vertexAt(line.numPoints()-1)
        return math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))

class GridSearchSolver:
    """
    Performs a brute-force grid search to find the optimal shift for maximal overlap.
    Range: -3.0m to +3.0m, Step: 0.25m
    """
    def __init__(self, current_geom, cadastral_geom, tolerance=0.1):
        self.current_geom = current_geom
        self.cadastral_geom = cadastral_geom
        self.tolerance = tolerance
        self.best_dx = 0.0
        self.best_dy = 0.0
        self.max_score = 0.0
        
    def solve(self):
        # Create Reference Buffer (Cadastral)
        # EndCapStyle: Flat for accurate line matching
        ref_buffer = self.cadastral_geom.buffer(self.tolerance, 8, Qgis.EndCapStyle.Flat, Qgis.JoinStyle.Round, 2.0)
        
        step = 0.25
        limit = 3.0
        
        current_len = self.current_geom.length()
        if current_len == 0: return 0.0, 0.0, 0.0
        
        best_score = -1.0
        best_dist_sq = float('inf') # Tie breaker: shortest shift
        
        # Ranges
        steps = int(limit / step)
        
        for x in range(-steps, steps + 1):
            dx = x * step
            for y in range(-steps, steps + 1):
                dy = y * step
                
                # Apply Shift
                shifted = QgsGeometry(self.current_geom)
                shifted.translate(dx, dy)
                
                # Check Overlap
                if not shifted.intersects(ref_buffer):
                    score = 0.0
                else:
                    intersection = shifted.intersection(ref_buffer)
                    score = intersection.length() / current_len
                
                # Update Best
                # Prioritize Score, then minimal shift
                if score > best_score:
                    best_score = score
                    self.best_dx = dx
                    self.best_dy = dy
                    self.max_score = score
                    best_dist_sq = dx*dx + dy*dy
                elif abs(score - best_score) < 0.001: # Check drift / Tie-break
                    dist_sq = dx*dx + dy*dy
                    if dist_sq < best_dist_sq:
                        self.best_dx = dx
                        self.best_dy = dy
                        best_dist_sq = dist_sq
                        
        return self.best_dx, self.best_dy, self.max_score

class ConfidenceStringMatcher:
    """
    Revised Matching Engine (REQ-2026-02-05).
    Supports Absolute Distance, Hybrid Buffer Check, and Grid Search.
    """
    def __init__(self, sigma=0.1):
        self.sigma = sigma

    def _get_bearing(self, p1, p2):
        """Calculates the bearing (0-180) between two points."""
        angle = math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))
        return angle % 180

    def _angle_diff(self, a1, a2):
        """Calculates the minimum difference between two bearings (0-180)."""
        diff = abs(a1 - a2) % 180
        return min(diff, 180 - diff)

    def process_pair(self, cur_feat, cad_feat, mode="simple", transform=None, tolerance=0.1):
        geom_curr = cur_feat.geometry()
        geom_cad = cad_feat.geometry()
        
        # Transform
        if transform:
            geom_curr.transform(transform)
            
        # 1. Absolute Distance (Hausdorff)
        h_dist = geom_curr.hausdorffDistance(geom_cad)
        
        result = {
            "status": "불부합",
            "score": h_dist,
            "nd_cost": 1.0,
            "shift_dx": 0.0,
            "shift_dy": 0.0,
            "rec_shift": "-",
            "reliability": "낮음"
        }
        
        # Pass Condition 1: Distance check
        if h_dist <= tolerance:
            result.update({
                "status": "부합",
                "rec_shift": "-",
                "reliability": "높음",
                "nd_cost": 0.0
            })
            return result
            
        # 2. Buffer Overlap Ratio (Hybrid check)
        cad_buffer = geom_cad.buffer(tolerance, 5)
        intersection_buf = geom_curr.intersection(cad_buffer)
        overlap_ratio = 0.0
        if geom_curr.length() > 0:
            overlap_ratio = intersection_buf.length() / geom_curr.length()
        
        result["nd_cost"] = 1.0 - overlap_ratio

        # 3. Interval-Based Sampling Analysis (REQ-2026-02-06)
        # Sample every 1 meter
        interval = 1.0
        densified_geom = geom_curr.densifyByDistance(interval)
        
        # Prepare comparison lines (Rings) from Cadastral geometry
        comparison_lines = []
        if geom_cad.type() == QgsWkbTypes.PolygonGeometry:
            if geom_cad.isMultipart():
                for poly in geom_cad.asMultiPolygon():
                    for ring in poly:
                        comparison_lines.append(QgsGeometry.fromPolylineXY(ring))
            else:
                for ring in geom_cad.asPolygon():
                    comparison_lines.append(QgsGeometry.fromPolylineXY(ring))
        else:
            comparison_lines.append(geom_cad)
            
        # Prepare comparison segments and their angles
        cad_segments = []
        for line_geom in comparison_lines:
            # Extract vertices to build segments (convert QgsPoint to QgsPointXY)
            pts = [QgsPointXY(v.x(), v.y()) for v in line_geom.vertices()]
            for j in range(len(pts) - 1):
                seg = QgsGeometry.fromPolylineXY([pts[j], pts[j+1]])
                angle = self._get_bearing(pts[j], pts[j+1])
                cad_segments.append((seg, angle))
            
        sample_vectors = []
        distances = []
        
        survey_vertices = [v for v in densified_geom.vertices()]
        for i in range(len(survey_vertices)):
            pt = QgsPointXY(survey_vertices[i].x(), survey_vertices[i].y())
            pt_geom = QgsGeometry.fromPointXY(pt)
            
            # 1. Determine local survey orientation
            if i < len(survey_vertices) - 1:
                s_angle = self._get_bearing(survey_vertices[i], survey_vertices[i+1])
            elif i > 0:
                s_angle = self._get_bearing(survey_vertices[i-1], survey_vertices[i])
            else:
                s_angle = 0 # Point geometry fallback (safety)

            # 2. Find best matching segment (nearest + aligned)
            best_pt = None
            min_d = float('inf')
            
            angle_tolerance = 45.0 # Increased tolerance to catch more 'roughly parallel' lines
            outlier_limit = 10.0 
            
            # Prioritized Search:
            # 1. Look for 'aligned' segments within outlier limits
            # 2. If none, look for 'any' segment within outlier limits (Fallback)
            
            candidate_aligned = None
            min_d_aligned = float('inf')
            
            candidate_any = None
            min_d_any = float('inf')
            
            for seg_geom, seg_angle in cad_segments:
                nr = seg_geom.nearestPoint(pt_geom)
                if not nr.isEmpty():
                    d = pt_geom.distance(nr)
                    
                    # Track absolute nearest (fallback candidate)
                    if d < min_d_any:
                        min_d_any = d
                        candidate_any = nr.asPoint()
                    
                    # Track aligned nearest (primary candidate)
                    if self._angle_diff(s_angle, seg_angle) < angle_tolerance:
                        if d < min_d_aligned:
                            min_d_aligned = d
                            candidate_aligned = nr.asPoint()

            # Decision Logic:
            # Prefer aligned candidate if it exists and is within rational distance
            # Otherwise use the absolute nearest candidate (e.g. at corners)
            
            if candidate_aligned and min_d_aligned < outlier_limit:
                best_pt = candidate_aligned
                min_d = min_d_aligned
            elif candidate_any and min_d_any < outlier_limit:
                 # Fallback: Angles didn't match, but this point is close enough to be considered
                best_pt = candidate_any
                min_d = min_d_any

            # Only record if it's within the outlier limit to prevent 'jumping'
            if best_pt and min_d < outlier_limit:
                distances.append(min_d)
                # IMPORTANT: Always append together to keep indices sync
                vec = QgsGeometry.fromPolylineXY([pt, best_pt])
                sample_vectors.append(vec)
        
        if distances:
            avg_dist = sum(distances) / len(distances)
            max_dist = max(distances)
        else:
            avg_dist = h_dist
            max_dist = h_dist

        if sample_vectors and len(sample_vectors) == len(distances):
            # Store all vectors as a MultiLineString for visualization
            result["error_vectors"] = QgsGeometry.fromMultiPolylineXY([v.asPolyline() for v in sample_vectors])
            # Correctly map max distance to its specific vector
            if distances:
                max_idx = distances.index(max(distances))
                if max_idx < len(sample_vectors):
                    result["error_line"] = sample_vectors[max_idx]
        
        result["score"] = avg_dist
        
        return result
        
    def process_multi_context(self, cur_feat, cad_feats, transform=None, tolerance=0.1):
        """
        1:N Context-Aware Matching (REQ-2026-02-06-TaskOrder).
        Matches current feature against multiple cadastral candidates using strict angle verification.
        """
        geom_curr = cur_feat.geometry()
        
        # Transform current geometry if needed
        if transform:
            geom_curr.transform(transform)

        # 1. Preprocess: Decompose ALL candidates into segments with angles
        target_segments = []
        for cf in cad_feats:
            c_geom = cf.geometry()
            # Handle Multi-parts and Polygons
            parts = []
            if c_geom.isMultipart():
                if c_geom.type() == QgsWkbTypes.PolygonGeometry:
                    for poly in c_geom.asMultiPolygon():
                        parts.extend([QgsGeometry.fromPolylineXY(ring) for ring in poly])
                else:
                    parts = c_geom.asMultiPolyline() # returns list of list of points? No, usually not geometry objects.
                    # PyQGIS asMultiPolyline returns [[pt, pt], [pt, pt]]
                    # We need QgsGeometry for easy vertex iteration or just iterate points
                    # Let's standardize on QgsGeometry for consistency with previous logic
                    temp_parts = c_geom.asMultiPolyline()
                    parts.extend([QgsGeometry.fromPolylineXY(p) for p in temp_parts])
            else:
                if c_geom.type() == QgsWkbTypes.PolygonGeometry:
                    for ring in c_geom.asPolygon():
                        parts.append(QgsGeometry.fromPolylineXY(ring))
                else:
                    parts.append(c_geom)

            # Build Segments
            for part_geom in parts:
                pts = [QgsPointXY(v.x(), v.y()) for v in part_geom.vertices()]
                for j in range(len(pts) - 1):
                    seg = QgsGeometry.fromPolylineXY([pts[j], pts[j+1]])
                    angle = self._get_bearing(pts[j], pts[j+1])
                    target_segments.append({'geom': seg, 'angle': angle})

        # 2. Interval Sampling & Matching
        interval = 1.0
        densified = geom_curr.densifyByDistance(interval)
        survey_pts = [QgsPointXY(v.x(), v.y()) for v in densified.vertices()]
        
        distances = []
        sample_vectors = []
        outlier_limit = 10.0
        angle_threshold = 30.0 # Strict parallel check

        for i, pt in enumerate(survey_pts):
            pt_geom = QgsGeometry.fromPointXY(pt)
            
            # Calculate Local Survey Angle
            # Use centered difference if possible, else forward/backward
            if i < len(survey_pts) - 1:
                s_angle = self._get_bearing(pt, survey_pts[i+1])
            elif i > 0:
                s_angle = self._get_bearing(survey_pts[i-1], pt)
            else:
                s_angle = 0
            
            # 3. Find Best Match among FILTERED segments
            best_pt = None
            min_d = float('inf')
            
            # Filter first by ANGLE
            candidates = []
            for seg in target_segments:
                diff = self._angle_diff(s_angle, seg['angle'])
                if diff <= angle_threshold:
                    candidates.append(seg['geom'])
            
            # Search nearest among candidates
            for cand_geom in candidates:
                nr = cand_geom.nearestPoint(pt_geom)
                if not nr.isEmpty():
                    d = nr.distance(pt_geom)
                    if d < min_d:
                        min_d = d
                        best_pt = nr.asPoint()
            
            # Valid Match Found?
            if best_pt and min_d < outlier_limit:
                distances.append(min_d)
                sample_vectors.append(QgsGeometry.fromPolylineXY([pt, best_pt]))
        
        # 4. Result Formatting
        if not distances:
            return {
                "status": "분석불가", "score": 999, "nd_cost": 1.0, 
                "reliability": "없음", "note": "매칭된 선분 없음"
            }

        avg_dist = sum(distances) / len(distances)
        max_dist = max(distances)
        
        result = {
            "status": "판단중",
            "score": avg_dist,
            "nd_cost": 0.0 if avg_dist < tolerance else 1.0, # Simplistic cost
            "reliability": "보통",
            "rec_shift": "-" 
        }

        # Status Logic
        if max_dist <= tolerance:
            result["status"] = "부합"
            result["reliability"] = "높음"
        else:
            result["status"] = "불부합"
        
        # Vectors
        if sample_vectors:
            result["error_vectors"] = QgsGeometry.fromMultiPolylineXY([v.asPolyline() for v in sample_vectors])
            max_idx = distances.index(max_dist)
            result["error_line"] = sample_vectors[max_idx]

        return result


    def calculate_parallel_overlap(self, geom1, geom2, angle_tolerance=10, dist_tolerance=2.0):
        """
        Calculates the total length of overlapping parallel segments between two lines.
        (REQ-2026-02-06-Hausdorff)
        """
        if geom1.type() != QgsWkbTypes.LineGeometry or geom2.type() != QgsWkbTypes.LineGeometry:
            return 0.0
            
        # Decompose geom1 into segments
        parts = []
        if geom1.isMultipart():
            parts = geom1.asMultiPolyline()
        else:
            parts = [geom1.asPolyline()]
            
        # Simplified approach: Densify geom1 and check closeness + parallelism to geom2.
        densified = geom1.densifyByDistance(0.5)
        d_pts = [QgsPointXY(p) for p in densified.vertices()]
        
        overlap_count = 0
        
        for i in range(len(d_pts)-1):
            p1 = d_pts[i]
            p2 = d_pts[i+1]
            mid_pt = QgsGeometry.fromPointXY(QgsPointXY((p1.x()+p2.x())/2, (p1.y()+p2.y())/2))
            
            # Check distance
            nearest = geom2.nearestPoint(mid_pt)
            if mid_pt.distance(nearest) > dist_tolerance:
                continue
                
            overlap_count += 1
            
        # Approx length
        return overlap_count * 0.5

class PointToLineAuditor:
    """
    P2L (Point-to-Line) Algorithm:
    1. Vertex Matching: Corners (>=45 deg) match to Cadastral Vertices.
    2. Contextual Matching: Short segments (<2m) inherit target from neighbors.
    3. Densification: Check every 2m.
    """
    def __init__(self, cadastral_feat, current_feat, densify_distance=2.0):
        self.cadastral_geom = cadastral_feat.geometry()
        self.current_geom = current_feat.geometry()
        self.densify_dist = densify_distance

    def process(self):
        # 1. Prepare Cadastral Data (Segments & Vertices)
        cad_segments_data = [] # List of (geometry, angle)
        cad_vertices = []
        
        # Robust extraction of lines from Polygon/Line geometries
        cad_parts = []
        if self.cadastral_geom.type() == QgsWkbTypes.PolygonGeometry:
            if self.cadastral_geom.isMultipart():
                for poly in self.cadastral_geom.asMultiPolygon():
                    cad_parts.extend(poly)
            else:
                cad_parts.extend(self.cadastral_geom.asPolygon())
        elif self.cadastral_geom.type() == QgsWkbTypes.LineGeometry:
            if self.cadastral_geom.isMultipart():
                cad_parts = self.cadastral_geom.asMultiPolyline()
            else:
                cad_parts = [self.cadastral_geom.asPolyline()]
        
        # Step 1: Reference Densification (for Vertices)
        # Apply densifyByDistance using the same distance as survey line
        densified_cad = self.cadastral_geom.densifyByDistance(self.densify_dist)
        for v in densified_cad.vertices():
            cad_vertices.append(QgsPointXY(v.x(), v.y()))
            
        for part in cad_parts:
            for i in range(len(part)):
                if i < len(part) - 1:
                    p1, p2 = part[i], part[i+1]
                    geom = QgsGeometry.fromPolylineXY([p1, p2])
                    # Calculate Azimuth (Angle)
                    angle = math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))
                    cad_segments_data.append((geom, angle))

        # 2. Process Current Geometry
        if self.current_geom.isMultipart():
            cur_parts = self.current_geom.asMultiPolyline()
        else:
            cur_parts = [self.current_geom.asPolyline()]

        max_dev = 0.0
        sum_sq = 0.0
        count = 0
        vectors = []
        sum_dx = 0.0
        sum_dy = 0.0
        
        for part in cur_parts:
            if len(part) < 2: continue
            
            # A. Identify Corners
            corners = [False] * len(part)
            for i in range(1, len(part) - 1):
                p_prev, p_curr, p_next = part[i-1], part[i], part[i+1]
                a1 = math.atan2(p_curr.y() - p_prev.y(), p_curr.x() - p_prev.x())
                a2 = math.atan2(p_next.y() - p_curr.y(), p_next.x() - p_curr.x())
                diff = math.degrees(abs(a1 - a2))
                if diff > 180: diff = 360 - diff
                if diff >= 45:
                    corners[i] = True
            
            # C. Iterative Check (Densified)
            for i in range(len(part) - 1):
                p_start, p_end = part[i], part[i+1]
                
                # Dynamic Target Selection: Find ALL parallel cadastral segments
                # This handles cases where the survey segment spans multiple cadastral segments
                curr_angle = math.degrees(math.atan2(p_end.y() - p_start.y(), p_end.x() - p_start.x()))
                
                candidates = []
                for geom, ang in cad_segments_data:
                    diff = abs(curr_angle - ang)
                    if diff > 180: diff = 360 - diff
                    if diff <= 30 or abs(diff - 180) <= 30: # Parallel check
                        candidates.append(geom)
                
                if candidates:
                    # Combine all parallel segments into one target geometry
                    target_geom = QgsGeometry.fromMultiPolylineXY([c.asPolyline() for c in candidates])
                else:
                    # Fallback to full cadastral geometry if no parallel found
                    target_geom = self.cadastral_geom
                
                seg_geom = QgsGeometry.fromPolylineXY([p_start, p_end])
                densified = seg_geom.densifyByDistance(self.densify_dist).asPolyline()
                
                # Process points (avoid double counting p_end unless last segment)
                points_to_check = densified if i == len(part) - 2 else densified[:-1]
                
                for j, pt in enumerate(points_to_check):
                    pt_geom = QgsGeometry.fromPointXY(pt)
                    is_corner = (j == 0 and corners[i])
                    if i == len(part) - 2 and j == len(points_to_check) - 1 and corners[i+1]:
                        is_corner = True
                        
                    if is_corner:
                        # Vertex Matching
                        near_pt = min(cad_vertices, key=lambda v: pt.sqrDist(v))
                        dist = math.sqrt(pt.sqrDist(near_pt))
                    else:
                        # Contextual Matching
                        # Step 3: Strict Perpendicular Projection
                        nearest_res = target_geom.nearestPoint(pt_geom)
                        near_pt = nearest_res.asPoint()
                        dist = pt_geom.distance(nearest_res)
            
                    if dist > max_dev: max_dev = dist
                    sum_sq += dist * dist
                    count += 1
                    sum_dx += (near_pt.x() - pt.x())
                    sum_dy += (near_pt.y() - pt.y())
                    
                    if dist >= 0.1:
                        vectors.append(QgsGeometry.fromPolylineXY([pt, near_pt]))
                
        rmse = math.sqrt(sum_sq / count) if count > 0 else 0.0
        avg_dx = sum_dx / count if count > 0 else 0.0
        avg_dy = sum_dy / count if count > 0 else 0.0
        
        multi_line = QgsGeometry.fromMultiPolylineXY([v.asPolyline() for v in vectors]) if vectors else QgsGeometry()
        
        return {
            "max_deviation": max_dev,
            "rmse": rmse,
            "avg_dx": avg_dx,
            "avg_dy": avg_dy,
            "error_vectors": multi_line
        }

class TopologyAuditor:
    """
    Advanced Topology Normalization:
    1. Split Process (1:N): Split survey lines by cadastral boundaries.
    2. Group Process (N:1): Group segments referencing the same cadastral edge.
    3. Weighted Average Error: Calculate error for the group.
    """
    def __init__(self, cad_layer, cur_layer, transform=None):
        self.cad_layer = cad_layer
        self.cur_layer = cur_layer
        self.xform = transform

    def process(self):
        results = []
        # 1. Build Spatial Index for Cadastral Layer
        cad_index = QgsSpatialIndex(self.cad_layer.getFeatures())
        
        # Intermediate storage for grouping
        # Key: target_edge_id (str), Value: dict of stats
        groups = {}
        
        # 2. Iterate & Split
        for cur_feat in self.cur_layer.getFeatures():
            geom = cur_feat.geometry()
            if not geom: continue
            if self.xform:
                geom.transform(self.xform)
                
            bbox = geom.boundingBox()
            cids = cad_index.intersects(bbox)
            
            for cid in cids:
                cad_feat = self.cad_layer.getFeature(cid)
                cad_geom = cad_feat.geometry()
                
                if not cad_geom.intersects(geom): continue
                
                intersection = geom.intersection(cad_geom)
                if intersection.isEmpty(): continue
                
                parts = intersection.asMultiPolyline() if intersection.isMultipart() else [intersection.asPolyline()]
                
                for part in parts:
                    if len(part) < 2: continue
                    seg_geom = QgsGeometry.fromPolylineXY(part)
                    
                    # Find Target Edge
                    target = self.find_target_edge(seg_geom, cad_feat)
                    if target:
                        tid = target['id']
                        if tid not in groups:
                            groups[tid] = {
                                'pnu': self.get_pnu(cad_feat),
                                'w_error_sum': 0.0,
                                'len_sum': 0.0,
                                'vectors': [],
                                'fids': []
                            }
                        
                        # Calculate average distance for this segment
                        d1 = target['geom'].distance(QgsGeometry.fromPointXY(QgsPointXY(part[0])))
                        d2 = target['geom'].distance(QgsGeometry.fromPointXY(QgsPointXY(part[-1])))
                        avg_dist = (d1 + d2) / 2.0
                        
                        groups[tid]['w_error_sum'] += avg_dist * seg_geom.length()
                        groups[tid]['len_sum'] += seg_geom.length()
                        groups[tid]['fids'].append(cur_feat.id())
                        
                        # Create vector for visualization (midpoint to line)
                        mid = seg_geom.interpolate(seg_geom.length()/2)
                        proj = target['geom'].nearestPoint(mid)
                        vec = QgsGeometry.fromPolylineXY([mid.asPoint(), proj.asPoint()])
                        groups[tid]['vectors'].append(vec)

        # 3. Finalize Groups
        for tid, g in groups.items():
            if g['len_sum'] == 0: continue
            avg_error = g['w_error_sum'] / g['len_sum']
            
            combined_vector = QgsGeometry()
            if g['vectors']:
                combined_vector = QgsGeometry.fromMultiPolylineXY([v.asPolyline() for v in g['vectors']])
                
            results.append({
                'pnu': g['pnu'],
                'error': avg_error,
                'length': g['len_sum'],
                'vector': combined_vector,
                'fid': g['fids'][0] if g['fids'] else None
            })
            
        return results

    def get_pnu(self, feat):
        for field in ['jibun', 'pnu', 'JIBUN', 'PNU']:
            if feat.fieldNameIndex(field) != -1:
                return str(feat[field])
        return str(feat.id())

    def find_target_edge(self, seg_geom, cad_feat):
        edges = []
        geom = cad_feat.geometry()
        if geom.isMultipart():
            polys = geom.asMultiPolygon()
        else:
            polys = [geom.asPolygon()]
            
        idx = 0
        for poly in polys:
            for ring in poly:
                for i in range(len(ring)-1):
                    p1, p2 = ring[i], ring[i+1]
                    edge = QgsGeometry.fromPolylineXY([p1, p2])
                    edges.append((f"{cad_feat.id()}_{idx}", edge))
                    idx += 1
        
        seg_angle = self._get_angle(seg_geom)
        best = None
        min_dist = float('inf')
        
        for eid, edge in edges:
            edge_angle = self._get_angle(edge)
            diff = abs(seg_angle - edge_angle)
            if diff > 180: diff = 360 - diff
            
            if diff <= 30 or abs(diff - 180) <= 30:
                dist = edge.distance(seg_geom)
                if dist < min_dist:
                    min_dist = dist
                    best = {'id': eid, 'geom': edge}
        
        return best

    def _get_angle(self, line):
        if line.isEmpty(): return 0
        p1 = line.vertexAt(0)
        p2 = line.vertexAt(line.numPoints()-1)
        return math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))