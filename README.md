# Loop over the input files
root -b script_LoopGenerators.C
root -b script_LoopGenerators_Weights.C

# To overlay the outputs of the loop w/o any Ac multiplication
# Generators and GENIE variations
root -b script_GeneratorOverlay.cxx
root -b GeneratorOverlay2D3D.cxx

# Interaction breakdown
root -b GeneratorInteBreakDown.cxx

# FSI Study
root -b GeneratorFSIStudy.cxx

# Graph overlay
root -b script_GraphOverlay.cxx

root -b script_PRD_3DXSec.C
root -b script_RatioGeneratorOverlay.cxx

# XSecs with Ac multiplication
root -b script_WienerSVD_OverlayGenerators.C
root -b script_TwoDimWienerSVD_OverlayGenerators.C
root -b script_Log_TwoDimWienerSVD_OverlayGenerators.C

