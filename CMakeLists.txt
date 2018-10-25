cmake_minimum_required( VERSION 3.2 FATAL_ERROR )

atlas_subdir( loop_the_loop )

atlas_depends_on_subdirs(
    PRIVATE
    Control/AthToolSupport/AsgTools
    Control/xAODRootAccess
    Event/xAOD/xAODCore
    Event/xAOD/xAODMuon
    Event/xAOD/xAODEgamma
    PhysicsAnalysis/Interfaces/EgammaAnalysisInterfaces
    PhysicsAnalysis/Interfaces/MuonAnalysisInterfaces
    PhysicsAnalysis/Interfaces/AsgAnalysisInterfaces
)

find_package( ROOT REQUIRED COMPONENTS RIO Hist Tree Net Core )
set(_common
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES ${ROOT_LIBRARIES}
    xAODRootAccess
    xAODMuon
    xAODEgamma
    AsgTools
    EgammaAnalysisInterfacesLib
)

atlas_add_executable( run_loop util/looper.cxx ${_common})
    
    
