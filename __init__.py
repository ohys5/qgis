# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Cadastral Auditor (지적 불부합 분석기)
                                 A QGIS plugin
 이 파일은 QGIS 플러그인의 진입점입니다.
 QGIS가 플러그인을 로드할 때 가장 먼저 이 파일을 찾습니다.
 ***************************************************************************/
"""

def classFactory(iface):
    """
    CadastralAuditor 클래스를 로드하고 인스턴스를 반환합니다.
    이 함수는 QGIS가 플러그인을 활성화할 때 호출됩니다.
    
    :param iface: QGIS 인터페이스에 대한 참조 (주 메뉴, 툴바 등 접근 가능)
    :return: CadastralAuditor 플러그인 인스턴스
    """
    # main_plugin.py 파일에서 CadastralAuditor 클래스를 가져옵니다.
    from .main_plugin import CadastralAuditor
    
    # 플러그인 인스턴스를 생성하여 QGIS에 반환합니다.
    return CadastralAuditor(iface)