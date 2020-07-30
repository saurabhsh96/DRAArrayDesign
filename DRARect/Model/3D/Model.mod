'# MWS Version: Version 2018.0 - Oct 26 2017 - ACIS 27.0.2 -

'# length = mm
'# frequency = GHz
'# time = ns
'# frequency range: fmin = 0.8*frequency_centre fmax = 1.2*frequency_centre


'@ use template: Antenna (Horn, Waveguide)

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
' Template for Antenna in Free Space
' ==================================
' (CSTxMWSxONLY)
' draw the bounding box
Plot.DrawBox True
' set units to mm, ghz
With Units 
     .Geometry "mm"
     .Frequency "ghz"
     .Time "ns" 
End With 
' set background material to vacuum
With Background 
     .Type "Normal" 
     .Epsilon "1.0" 
     .Mue "1.0" 
     .XminSpace "0.0" 
     .XmaxSpace "0.0" 
     .YminSpace "0.0" 
     .YmaxSpace "0.0" 
     .ZminSpace "0.0" 
     .ZmaxSpace "0.0" 
End With 
' set boundary conditions to open
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "none" 
     .Ysymmetry "none" 
     .Zsymmetry "none" 
End With
' switch on FD-TET setting for accurate farfields
FDSolver.ExtrudeOpenBC "True" 
Mesh.FPBAAvoidNonRegUnite "True" 
Mesh.ConsiderSpaceForLowerMeshLimit "False" 
Mesh.MinimumStepNumber "5"

'@ new component: antenna

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Component.New "antenna"

'@ define brick: antenna:groundplane

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Brick
     .Reset 
     .Name "groundplane" 
     .Component "antenna" 
     .Material "PEC" 
     .Xrange "-dielectric_width*2", "dielectric_width*2" 
     .Yrange "-dielectric_length*2", "dielectric_length*2" 
     .Zrange "-metal_thickness", "0" 
     .Create
End With

'@ define material: dielectric

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Material 
     .Reset 
     .Name "dielectric"
     .FrqType "all" 
     .Type "Normal" 
     .Epsilon "relative_permittivity" 
     .Mue "1" 
     .Kappa "0" 
     .TanD "0.0" 
     .TanDFreq "0.0" 
     .TanDGiven "False" 
     .TanDModel "ConstTanD" 
     .KappaM "0" 
     .TanDM "0.0" 
     .TanDMFreq "0.0" 
     .TanDMGiven "False" 
     .TanDMModel "ConstTanD" 
     .DispModelEps "None" 
     .DispModelMue "None" 
     .DispersiveFittingSchemeEps "General 1st" 
     .DispersiveFittingSchemeMue "General 1st" 
     .UseGeneralDispersionEps "False" 
     .UseGeneralDispersionMue "False" 
     .Rho "0" 
     .ThermalType "Normal" 
     .ThermalConductivity "0" 
     .HeatCapacity "0" 
     .MetabolicRate "0" 
     .BloodFlow "0" 
     .Colour "0", "1", "1" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .Create
End With

'@ define brick: antenna:DR

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Brick
     .Reset 
     .Name "DR" 
     .Component "antenna" 
     .Material "dielectric" 
     .Xrange "-dielectric_width/2", "dielectric_width/2" 
     .Yrange "-dielectric_length/2", "dielectric_length/2" 
     .Zrange "0", "dielectric_height" 
     .Create
End With

'@ define cylinder: antenna:probe_cutout

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Cylinder 
     .Reset 
     .Name "probe_cutout" 
     .Component "antenna" 
     .Material "PEC" 
     .OuterRadius "coax_diameter/2" 
     .InnerRadius "0" 
     .Axis "z" 
     .Zrange "-metal_thickness", "0" 
     .Xcenter "dielectric_width/2-probe_inset" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ boolean subtract shapes: antenna:groundplane, antenna:probe_cutout

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Solid 
     .Version 9
     .Subtract "antenna:groundplane", "antenna:probe_cutout" 
End With

'@ define cylinder: antenna:probe

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Cylinder 
     .Reset 
     .Name "probe" 
     .Component "antenna" 
     .Material "PEC" 
     .OuterRadius "probe_diameter/2" 
     .InnerRadius "0" 
     .Axis "z" 
     .Zrange "-metal_thickness-coax_length", "probe_height" 
     .Xcenter "dielectric_width/2-probe_inset" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ define cylinder: antenna:coax

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Cylinder 
     .Reset 
     .Name "coax" 
     .Component "antenna" 
     .Material "PEC" 
     .OuterRadius "coax_diameter/2+metal_thickness" 
     .InnerRadius "coax_diameter/2" 
     .Axis "z" 
     .Zrange "-metal_thickness-coax_length", "-metal_thickness" 
     .Xcenter "dielectric_width/2-probe_inset" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ define frequency range

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Solver.FrequencyRange "0.8*frequency_centre", "1.2*frequency_centre"

'@ define boundaries

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "none" 
     .Ysymmetry "magnetic" 
     .Zsymmetry "none" 
     .XminThermal "isothermal" 
     .XmaxThermal "isothermal" 
     .YminThermal "isothermal" 
     .YmaxThermal "isothermal" 
     .ZminThermal "isothermal" 
     .ZmaxThermal "isothermal" 
     .XsymmetryThermal "none" 
     .YsymmetryThermal "isothermal" 
     .ZsymmetryThermal "none" 
     .ApplyInAllDirections "False" 
     .XminTemperature "" 
     .XminTemperatureType "None" 
     .XmaxTemperature "" 
     .XmaxTemperatureType "None" 
     .YminTemperature "" 
     .YminTemperatureType "None" 
     .YmaxTemperature "" 
     .YmaxTemperatureType "None" 
     .ZminTemperature "" 
     .ZminTemperatureType "None" 
     .ZmaxTemperature "" 
     .ZmaxTemperatureType "None" 
End With

'@ pick edge

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Pick.PickEdgeFromId "antenna:coax", "2", "2"

'@ define port: 1

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Port 
     .Reset 
     .PortNumber "1" 
     .NumberOfModes "1" 
     .AdjustPolarization False 
     .PolarizationAngle "0.0" 
     .ReferencePlaneDistance "0.0" 
     .TextSize "50" 
     .Coordinates "Picks" 
     .Orientation "positive" 
     .PortOnBound "False" 
     .ClipPickedPortToBound "False" 
     .Xrange "0.50667123963133", "0.67715276036867" 
     .Yrange "-0.085240760368665", "0.085240760368665" 
     .Zrange "-0.2369565624", "-0.2369565624" 
     .XrangeAdd "0.0", "0.0" 
     .YrangeAdd "0.0", "0.0" 
     .ZrangeAdd "0.0", "0.0" 
     .SingleEnded "False" 
     .Create 
End With

'@ pick face

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Pick.PickFaceFromId "port1", "1"

'@ define material: PEC_clear

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Material 
     .Reset 
     .Name "PEC_clear"
     .FrqType "all" 
     .Type "Pec" 
     .Epsilon "1.0" 
     .Mue "1.0" 
     .Rho "0" 
     .ThermalType "Normal" 
     .ThermalConductivity "0" 
     .HeatCapacity "0" 
     .MetabolicRate "0" 
     .BloodFlow "0" 
     .Colour "0.647059", "0.666667", "0.72549" 
     .Wireframe "True" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .Create
End With

'@ define extrude: antenna:port_back

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Extrude 
     .Reset 
     .Name "port_back" 
     .Component "antenna" 
     .Material "PEC_clear" 
     .Mode "Picks" 
     .Height "-coax_length/4" 
     .Twist "0.0" 
     .Taper "0.0" 
     .UsePicksForHeight "False" 
     .DeleteBaseFaceSolid "False" 
     .ClearPickedFace "True" 
     .Create 
End With

'@ define automesh for: antenna:probe

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Solid.SetSolidMeshProperties "antenna:probe", "PBA", "True", "0", "True", "True", "probe_diameter/4", "probe_diameter/4", "0", "probe_diameter*2", "probe_diameter*2", "0", "False", "1.0", "False", "1.0", "False", "True", "True"

'@ define farfield monitor: farfield (f=frequency_centre)

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Monitor 
     .Reset 
     .Name "farfield (f=frequency_centre)" 
     .Domain "Frequency" 
     .FieldType "Farfield" 
	    .Frequency "frequency_centre" 
     .Create 
End With

'@ define solver parameters

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .CalculationType "TD-S" 
     .StimulationPort "All" 
     .StimulationMode "All" 
     .SteadyStateLimit "-40" 
     .MeshAdaption "True" 
     .AutoNormImpedance "True" 
     .NormingImpedance "input_resistance" 
     .CalculateModesOnly "False" 
     .SParaSymmetry "False" 
     .StoreTDResultsInCache "False" 
     .FullDeembedding "False" 
     .UseDistributedComputing "False" 
     .DistributeMatrixCalculation "False" 
     .MPIParallelization "False" 
     .SuperimposePLWExcitation "False" 
End With

'@ set 3d mesh adaptation results

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Mesh 
    .LinesPerWavelength "25" 
    .MinimumStepNumber "20" 
End With

'@ deactivate transient solver mesh adaptation

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Solver.MeshAdaption "False"

'@ define solver parameters

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .CalculationType "TD-S" 
     .StimulationPort "All" 
     .StimulationMode "All" 
     .SteadyStateLimit "-40" 
     .MeshAdaption "False" 
     .AutoNormImpedance "True" 
     .NormingImpedance "input_resistance" 
     .CalculateModesOnly "False" 
     .SParaSymmetry "False" 
     .StoreTDResultsInCache "False" 
     .FullDeembedding "False" 
     .UseDistributedComputing "False" 
     .DistributeMatrixCalculation "False" 
     .MPIParallelization "False" 
     .SuperimposePLWExcitation "False" 
End With

'@ define automesh for: antenna:probe

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
Solid.SetSolidMeshProperties "antenna:probe", "PBA", "True", "0", "True", "True", "probe_diameter/6", "probe_diameter/6", "0", "probe_inset-probe_diameter/2", "probe_diameter*4", "0", "False", "1.0", "False", "1.0", "False", "True", "True"

'@ set mesh properties

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Mesh 
     .UseRatioLimit "True" 
     .RatioLimit "20" 
     .LinesPerWavelength "20" 
     .MinimumStepNumber "20" 
     .Automesh "True" 
     .MeshType "PBA" 
     .SetCreator "High Frequency" 
End With

'@ define boundaries

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "none" 
     .Ysymmetry "none" 
     .Zsymmetry "none" 
     .XminThermal "isothermal" 
     .XmaxThermal "isothermal" 
     .YminThermal "isothermal" 
     .YmaxThermal "isothermal" 
     .ZminThermal "isothermal" 
     .ZmaxThermal "isothermal" 
     .XsymmetryThermal "none" 
     .YsymmetryThermal "none" 
     .ZsymmetryThermal "none" 
     .ApplyInAllDirections "False" 
     .XminTemperature "" 
     .XminTemperatureType "None" 
     .XmaxTemperature "" 
     .XmaxTemperatureType "None" 
     .YminTemperature "" 
     .YminTemperatureType "None" 
     .YmaxTemperature "" 
     .YmaxTemperatureType "None" 
     .ZminTemperature "" 
     .ZminTemperatureType "None" 
     .ZmaxTemperature "" 
     .ZmaxTemperatureType "None" 
End With

'@ set mesh properties

'[VERSION]2009.7|18.0.3|20090230[/VERSION]
With Mesh 
     .UseRatioLimit "True" 
     .RatioLimit "40" 
     .LinesPerWavelength "40" 
     .MinimumStepNumber "10" 
     .Automesh "True" 
     .MeshType "PBA" 
     .SetCreator "High Frequency" 
End With

'@ set units in materials

'[VERSION]2010.5|20.0.0|20100711[/VERSION]
Material.SetUnitInMaterial "PEC_clear", "GHz", "mm" 
Material.SetUnitInMaterial "dielectric", "GHz", "mm"

'@ switch working plane

'[VERSION]2010.5|20.0.0|20100711[/VERSION]
Plot.DrawWorkplane "false"

'@ define solver parameters

'[VERSION]2010.5|20.0.0|20100711[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .CalculationType "TD-S"
     .StimulationPort "All"
     .StimulationMode "All"
     .SteadyStateLimit "-40"
     .MeshAdaption "False"
     .AutoNormImpedance "True"
     .NormingImpedance "input_resistance"
     .CalculateModesOnly "False"
     .SParaSymmetry "False"
     .StoreTDResultsInCache  "False"
     .FullDeembedding "False"
     .SuperimposePLWExcitation "False"
End With

'@ define frequency domain solver parameters

'[VERSION]2010.5|20.0.0|20100711[/VERSION]
Mesh.SetCreator "High Frequency" 
With FDSolver
     .Reset 
     .Method "Tetrahedral Mesh" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "True" 
     .NormingImpedance "input_resistance" 
     .ModesOnly "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .CalculateExcitationsInParallel "True" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "False" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .HexMORSettings "", "1001" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "True" 
     .SetOpenBCTypeHex "Default" 
     .AddMonitorSamples "True" 
     .SParameterSweep "True" 
     .CalcStatBField "False" 
     .UseDoublePrecision "False" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "true" 
     .UseFastRCSSweepIntEq "true" 
     .UseSensitivityAnalysis "False" 
     .SetRCSSweepProperties "0.0", "0.0", "0","0.0", "0.0", "0", "0" 
     .SweepErrorThreshold "True", "0.01" 
     .SweepErrorChecks "2" 
     .SweepMinimumSamples "3" 
     .SweepConsiderAll "True" 
     .SweepConsiderReset 
     .InterpolationSamples "1001" 
     .SweepWeightEvanescent "1.0" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .LimitCPUs "True"
     .MaxCPUs "16"
End With
With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .PreconditionerType "Auto" 
End With
With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
End With

'@ set tetrahedral mesh type

'[VERSION]2010.5|20.0.0|20100711[/VERSION]
Mesh.MeshType "Tetrahedral"

'@ set units in materials

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
Material.SetUnitInMaterial "$CoilMaterial$", "GHz", "mm"

'@ set mesh properties (for backward compatibility)

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
With MeshSettings
     .SetMeshType "Hex"
     .Set "Version", 0%
     .SetMeshType "Tet"
     .Set "Version", 0%
     .SetMeshType "Srf"
     .Set "Version", 0%
End With
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "CellsPerWavelengthPolicy", "cellsperwavelength" 
     .Set "CurvatureOrderPolicy", "off" 
     .SetMeshType "Plane" 
     .Set "CurvatureOrderPolicy", "off" 
End With

'@ change solver and mesh type

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
ChangeSolverAndMeshType "HF Frequency Domain"

'@ define frequency domain solver parameters

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
Mesh.SetCreator "High Frequency" 
With FDSolver
     .Reset 
     .SetMethod "Tetrahedral", "General purpose" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "True" 
     .NormingImpedance "input_resistance" 
     .ModesOnly "False" 
     .ConsiderPortLossesTet "True" 
     .SetShieldAllPorts "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .SetCalculateExcitationsInParallel "True", "False", "" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "False" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .FreqDistAdaptMode "Distributed" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "True" 
     .SetOpenBCTypeHex "Default" 
     .SetOpenBCTypeTet "Default" 
     .AddMonitorSamples "True" 
     .CalcStatBField "False" 
     .CalcPowerLoss "True" 
     .CalcPowerLossPerComponent "False" 
     .StoreSolutionCoefficients "True" 
     .UseDoublePrecision "False" 
     .UseDoublePrecision_ML "True" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "true" 
     .UseFastRCSSweepIntEq "true" 
     .UseSensitivityAnalysis "False" 
     .SweepErrorThreshold "True", "0.01" 
     .SweepErrorChecks "2" 
     .SweepMinimumSamples "3" 
     .SweepConsiderAll "True" 
     .SweepConsiderReset 
     .SetNumberOfResultDataSamples "1001" 
     .SweepWeightEvanescent "1.0" 
     .AccuracyROM "1e-4" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .LimitCPUs "True"
     .MaxCPUs "16"
End With
With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .PreconditionerType "Auto" 
End With
With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
     .DeembedExternalPorts "True" 
     .SetOpenBC_XY "True" 
     .OldRCSSweepDefintion "True" 
     .SetAccuracySetting "Custom" 
End With

'@ set mesh properties (Tetrahedral)

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
With Mesh 
     .MeshType "Tetrahedral" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "4" 
     .Set "StepsPerWaveFar", "4" 
     .Set "PhaseErrorNear", "0.02" 
     .Set "PhaseErrorFar", "0.02" 
     .Set "CellsPerWavelengthPolicy", "cellsperwavelength" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "10" 
     .Set "StepsPerBoxFar", "1" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     'MIN CELL 
     .Set "UseRatioLimit", "0" 
     .Set "RatioLimit", "100" 
     .Set "MinStep", "0" 
     'MESHING METHOD 
     .SetMeshType "Unstr" 
     .Set "Method", "0" 
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "CurvatureOrder", "1" 
     .Set "CurvatureOrderPolicy", "off" 
     .Set "CurvRefinementControl", "NormalTolerance" 
     .Set "NormalTolerance", "22.5" 
     .Set "SrfMeshGradation", "2" 
     .Set "SrfMeshOptimization", "1" 
End With 
With MeshSettings 
     .SetMeshType "Unstr" 
     .Set "UseMaterials",  "1" 
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "UseAnisoCurveRefinement", "1" 
     .Set "UseSameSrfAndVolMeshGradation", "1" 
     .Set "VolMeshGradation", "2" 
     .Set "VolMeshOptimization", "1" 
End With 
With MeshSettings 
     .SetMeshType "Unstr" 
     .Set "SmallFeatureSize", "0" 
     .Set "CoincidenceTolerance", "1e-006" 
     .Set "SelfIntersectionCheck", "1" 
End With 
With MeshSettings 
     .SetMeshType "Unstr" 
     .Set "UseDC", "0" 
End With 
With Mesh 
     .SetParallelMesherMode "Tet", "maximum" 
     .SetMaxParallelMesherThreads "Tet", "1" 
End With

'@ delete shape: antenna:port_back

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
Solid.Delete "antenna:port_back"

'@ Exported from Antenna Magus: Rectangular dielectric resonator antenna (DRA) - Friday, April 10, 2020

'[VERSION]2014.6|23.0.0|20090230[/VERSION]

'@ create group: meshgroup1

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
Group.Add "meshgroup1", "mesh"


'@ set local mesh properties for: meshgroup1

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With MeshSettings
     With .ItemMeshSettings ("group$meshgroup1")
          .SetMeshType "Hex"
          .Set "EdgeRefinement", "1.0"
          .Set "Extend", "probe_inset-probe_diameter/2", "probe_diameter*4", "0"
          .Set "Fixpoints", 1
          .Set "MeshType", "Default"
          .Set "NumSteps", 0, 0, 0
          .Set "Priority", "0"
          .Set "RefinementPolicy", "ABS_VALUE"
          .Set "SnappingIntervals", 0, 0, 0
          .Set "SnappingPriority", 0
          .Set "SnapTo", 1, 1, 1
          .Set "Step", "probe_diameter/6", "probe_diameter/6", "0"
          .Set "StepRatio", 0, 0, 0
          .Set "StepRefinementCollectPolicy", "REFINE_ALL"
          .Set "StepRefinementExtentPolicy", "EXTENT_ABS_VALUE"
          .Set "UseDielectrics", 1
          .Set "UseEdgeRefinement", 0
          .Set "UseForRefinement", 1
          .Set "UseForSnapping", 1
          .Set "UseSameExtendXYZ", 1
          .Set "UseSameStepWidthXYZ", 1
          .Set "UseSnappingPriority", 0
          .Set "UseStepAndExtend", 1
          .Set "UseVolumeRefinement", 0
          .Set "VolumeRefinement", "1.0"
          .SetMeshType "HexTLM"
          .Set "EdgeRefinement", "1.0"
          .Set "Extend", "probe_inset-probe_diameter/2", "probe_diameter*4", "0"
          .Set "NumSteps", 0, 0, 0
          .Set "RefinementPolicy", "ABS_VALUE"
          .Set "SnappingIntervals", 0, 0, 0
          .Set "SnappingPriority", 0
          .Set "SnapTo", 1, 1, 1
          .Set "Step", "probe_diameter/6", "probe_diameter/6", "0"
          .Set "StepRatio", 0, 0, 0
          .Set "StepRefinementCollectPolicy", "REFINE_ALL"
          .Set "StepRefinementExtentPolicy", "EXTENT_ABS_VALUE"
          .Set "UnlumpWithin", 0
          .Set "UseDielectrics", 1
          .Set "UseEdgeRefinement", 0
          .Set "UseForRefinement", 1
          .Set "UseForSnapping", 1
          .Set "UseSameExtendXYZ", 1
          .Set "UseSameStepWidthXYZ", 1
          .Set "UseSnappingPriority", 0
          .Set "UseStepAndExtend", 1
          .Set "UseVolumeRefinement", 0
          .Set "VolumeRefinement", "1.0"
     End With
End With


'@ add items to group: "meshgroup1"

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
Group.AddItem "solid$antenna:probe", "meshgroup1"


'@ set local mesh properties (for backward compatibility)

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With MeshSettings
     With .ItemMeshSettings ("group$meshgroup1")
          .SetMeshType "Hex"
          .Set "UseSameStepWidthXYZ", 0
          .SetMeshType "Hex"
          .Set "UseSameExtendXYZ", 0
          .SetMeshType "HexTLM"
          .Set "UseSameStepWidthXYZ", 0
          .SetMeshType "HexTLM"
          .Set "UseSameExtendXYZ", 0
     End With
End With


'@ set 3d mesh adaptation properties (for backward compatibility)

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With MeshAdaption3D
    .SetType "HighFrequencyTet" 
    .SingularEdgeRefinement "0" 
    .MaxDeltaS "0.01" 
    .NumberOfDeltaSChecks "2" 
End With


'@ create legacy 1D signals (for backward compatibility)

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
FDSolver.CreateLegacy1DSignals "True" 


'@ define frequency domain solver parameters

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
Mesh.SetCreator "High Frequency" 

With FDSolver
     .Reset 
     .SetMethod "Tetrahedral", "General purpose" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "True" 
     .NormingImpedance "input_resistance" 
     .ModesOnly "False" 
     .ConsiderPortLossesTet "True" 
     .SetShieldAllPorts "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .SetCalcBlockExcitationsInParallel "True", "True", "" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "False" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .FreqDistAdaptMode "Distributed" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "True" 
     .SetOpenBCTypeHex "Default" 
     .SetOpenBCTypeTet "Default" 
     .AddMonitorSamples "True" 
     .CalcStatBField "False" 
     .CalcPowerLoss "True" 
     .CalcPowerLossPerComponent "False" 
     .StoreSolutionCoefficients "True" 
     .UseDoublePrecision "False" 
     .UseDoublePrecision_ML "True" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "True" 
     .UseFastRCSSweepIntEq "true" 
     .UseSensitivityAnalysis "False" 
     .RemoveAllStopCriteria "Hex"
     .AddStopCriterion "All S-Parameters", "0.01", "2", "Hex", "True"
     .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Hex", "False"
     .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Hex", "False"
     .RemoveAllStopCriteria "Tet"
     .AddStopCriterion "All S-Parameters", "0.01", "2", "Tet", "True"
     .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Tet", "False"
     .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Tet", "False"
     .AddStopCriterion "All Probes", "0.05", "2", "Tet", "True"
     .RemoveAllStopCriteria "Srf"
     .AddStopCriterion "All S-Parameters", "0.01", "2", "Srf", "True"
     .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Srf", "False"
     .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Srf", "False"
     .SweepMinimumSamples "3" 
     .SetNumberOfResultDataSamples "1001" 
     .SetResultDataSamplingMode "Automatic" 
     .SweepWeightEvanescent "1.0" 
     .AccuracyROM "1e-4" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .CreateLegacy1DSignals "True" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .UseParallelization "True"
     .MaxCPUs "16"
     .MaximumNumberOfCPUDevices "2"
End With

With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .SetRealGroundMaterialName "" 
     .CalcFarFieldInRealGround "False" 
     .RealGroundModelType "Auto" 
     .PreconditionerType "Auto" 
     .ExtendThinWireModelByWireNubs "False" 
End With

With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
     .DeembedExternalPorts "True" 
     .SetOpenBC_XY "True" 
     .OldRCSSweepDefintion "True" 
     .SetAccuracySetting "Custom" 
     .CalculateSParaforFieldsources "True" 
     .ModeTrackingCMA "True" 
     .NumberOfModesCMA "3" 
     .StartFrequencyCMA "-1.0" 
     .SetAccuracySettingCMA "Default" 
     .FrequencySamplesCMA "0" 
     .SetMemSettingCMA "Auto" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ change solver type

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
ChangeSolverType "HF Frequency Domain" 


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Efield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Efield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.692307692308", "-57.692307692308", "0", "0.00053499842393752", "0" 
     .Antenna "-34.615384615385", "-57.692307692308", "0", "0.0021046954042014", "0" 
     .Antenna "-11.538461538461", "-57.692307692308", "0", "0.0058898004886415", "0" 
     .Antenna "11.538461538461", "-57.692307692308", "0", "0.013559865528552", "0" 
     .Antenna "34.615384615385", "-57.692307692308", "0", "0.027350796859432", "0" 
     .Antenna "57.692307692308", "-57.692307692308", "0", "0.049934024036287", "0" 
     .Antenna "-57.692307692308", "-34.615384615385", "0", "0.084143679636849", "0" 
     .Antenna "-34.615384615385", "-34.615384615385", "0", "0.13255937006707", "0" 
     .Antenna "-11.538461538461", "-34.615384615385", "0", "0.19698705148557", "0" 
     .Antenna "11.538461538461", "-34.615384615385", "0", "0.27791621645135", "0" 
     .Antenna "34.615384615385", "-34.615384615385", "0", "0.37405671014422", "0" 
     .Antenna "57.692307692308", "-34.615384615385", "0", "0.48206447283712", "0" 
     .Antenna "-57.692307692308", "-11.538461538461", "0", "0.5965473658611", "0" 
     .Antenna "-34.615384615385", "-11.538461538461", "0", "0.71040047342638", "0" 
     .Antenna "-11.538461538461", "-11.538461538461", "0", "0.81546144893755", "0" 
     .Antenna "11.538461538461", "-11.538461538461", "0", "0.90341221853988", "0" 
     .Antenna "34.615384615385", "-11.538461538461", "0", "0.96679785244378", "0" 
     .Antenna "57.692307692308", "-11.538461538461", "0", "1", "0" 
     .Antenna "-57.692307692308", "11.538461538461", "0", "1", "0" 
     .Antenna "-34.615384615385", "11.538461538461", "0", "0.96679785244378", "0" 
     .Antenna "-11.538461538461", "11.538461538461", "0", "0.90341221853988", "0" 
     .Antenna "11.538461538461", "11.538461538461", "0", "0.81546144893755", "0" 
     .Antenna "34.615384615385", "11.538461538461", "0", "0.71040047342638", "0" 
     .Antenna "57.692307692308", "11.538461538461", "0", "0.5965473658611", "0" 
     .Antenna "-57.692307692308", "34.615384615385", "0", "0.48206447283712", "0" 
     .Antenna "-34.615384615385", "34.615384615385", "0", "0.37405671014422", "0" 
     .Antenna "-11.538461538461", "34.615384615385", "0", "0.27791621645135", "0" 
     .Antenna "11.538461538461", "34.615384615385", "0", "0.19698705148557", "0" 
     .Antenna "34.615384615385", "34.615384615385", "0", "0.13255937006707", "0" 
     .Antenna "57.692307692308", "34.615384615385", "0", "0.084143679636849", "0" 
     .Antenna "-57.692307692308", "57.692307692308", "0", "0.049934024036287", "0" 
     .Antenna "-34.615384615385", "57.692307692308", "0", "0.027350796859432", "0" 
     .Antenna "-11.538461538461", "57.692307692308", "0", "0.013559865528552", "0" 
     .Antenna "11.538461538461", "57.692307692308", "0", "0.0058898004886415", "0" 
     .Antenna "34.615384615385", "57.692307692308", "0", "0.0021046954042014", "0" 
     .Antenna "57.692307692308", "57.692307692308", "0", "0.00053499842393752", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.000535", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.0021047", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.0058898", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.01355987", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.0273508", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.04993402", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.08414368", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.13255937", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.19698705", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.27791622", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.37405671", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.48206447", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.59654737", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.71040047", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "0.81546145", "0" 
     .Antenna "11.53846", "-11.53846", "0", "0.90341222", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.96679785", "0" 
     .Antenna "57.69231", "-11.53846", "0", "1", "0" 
     .Antenna "-57.69231", "11.53846", "0", "1", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.96679785", "0" 
     .Antenna "-11.53846", "11.53846", "0", "0.90341222", "0" 
     .Antenna "11.53846", "11.53846", "0", "0.81546145", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.71040047", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.59654737", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.48206447", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.37405671", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.27791622", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.19698705", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.13255937", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.08414368", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.04993402", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.0273508", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.01355987", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.0058898", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.0021047", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.000535", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.15311473", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.16383204", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.17433372", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.18454487", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.19439152", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.20380132", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.21270429", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.22103349", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.22872566", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.23572192", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.24196827", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.24741623", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.25202327", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "0.25575325", "0" 
     .Antenna "11.53846", "-11.53846", "0", "0.25857683", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.26047171", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.26142292", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.26142292", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.26047171", "0" 
     .Antenna "-11.53846", "11.53846", "0", "0.25857683", "0" 
     .Antenna "11.53846", "11.53846", "0", "0.25575325", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.25202327", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.24741623", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.24196827", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.23572192", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.22872566", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.22103349", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.21270429", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.20380132", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.19439152", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.18454487", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.17433372", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.16383204", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.15311473", "0" 
     .Antenna "57.69231", "57.69231", "0", "1", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.09430874", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.09837463", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.10229254", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.10604379", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.10961036", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.11297499", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.11612126", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.11903376", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.1216981", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.1241011", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.12623078", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.12807651", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.12962904", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "0.13088055", "0" 
     .Antenna "11.53846", "-11.53846", "0", "0.13182476", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.1324569", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.13277377", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.13277377", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.1324569", "0" 
     .Antenna "-11.53846", "11.53846", "0", "0.13182476", "0" 
     .Antenna "11.53846", "11.53846", "0", "0.13088055", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.12962904", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.12807651", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.12623078", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.1241011", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.1216981", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.11903376", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.11612126", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.11297499", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.10961036", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.10604379", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.10229254", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.09837463", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.09430874", "0" 
     .Antenna "57.69231", "57.69231", "0", "1", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "8", "wavelength_centre/2", "0" 
     .YSet "8", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.60712017", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.68083915", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.68083915", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.60712017", "0" 
     .Antenna "57.69231", "-57.69231", "0", "1", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.60712017", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.3685949", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.41335118", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.41335118", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.3685949", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.60712017", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.68083915", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.41335118", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "0.46354194", "0" 
     .Antenna "11.53846", "-11.53846", "0", "0.46354194", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.41335118", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.68083915", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.68083915", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.41335118", "0" 
     .Antenna "-11.53846", "11.53846", "0", "0.46354194", "0" 
     .Antenna "11.53846", "11.53846", "0", "0.46354194", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.41335118", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.68083915", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.60712017", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.3685949", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.41335118", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.41335118", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.3685949", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.60712017", "0" 
     .Antenna "-57.69231", "57.69231", "0", "1", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.60712017", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.68083915", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.68083915", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.60712017", "0" 
     .Antenna "57.69231", "57.69231", "0", "1", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.40089934", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.50562471", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.63316612", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.63316612", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.50562471", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.40089934", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.50562471", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.63770708", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.79856564", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.79856564", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.63770708", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.50562471", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.63316612", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.79856564", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.79856564", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.63316612", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.63316612", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.79856564", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.79856564", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.63316612", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.50562471", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.63770708", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.79856564", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.79856564", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.63770708", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.50562471", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.40089934", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.50562471", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.63316612", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.63316612", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.50562471", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.40089934", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.14935406", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.28084414", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.38646353", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.38646353", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.28084414", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.14935406", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.28084414", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.52809698", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.72670281", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.72670281", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.52809698", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.28084414", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.38646353", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.72670281", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.72670281", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.38646353", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.38646353", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.72670281", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.72670281", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.38646353", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.28084414", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.52809698", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.72670281", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.72670281", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.52809698", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.28084414", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.14935406", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.28084414", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.38646353", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.38646353", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.28084414", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.14935406", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.02430759", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.08970269", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.15590893", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.15590893", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.08970269", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.02430759", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.08970269", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.3310312", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.57535311", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.57535311", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.3310312", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.08970269", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.15590893", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.57535311", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.57535311", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.15590893", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.15590893", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.57535311", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.57535311", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.15590893", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.08970269", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.3310312", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.57535311", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.57535311", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.3310312", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.08970269", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.02430759", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.08970269", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.15590893", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.15590893", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.08970269", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.02430759", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0175963", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.0726541", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.13265103", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.13265103", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.0726541", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0175963", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.0726541", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.29998465", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.54770854", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.54770854", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.29998465", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.0726541", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.13265103", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.54770854", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.54770854", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.13265103", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.13265103", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.54770854", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.54770854", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.13265103", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.0726541", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.29998465", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.54770854", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.54770854", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.29998465", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.0726541", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0175963", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.0726541", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.13265103", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.13265103", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.0726541", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0175963", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01568694", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06737416", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.12524751", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.12524751", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06737416", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.01568694", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06737416", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28936672", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53792817", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.53792817", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.28936672", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06737416", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.12524751", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53792817", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.53792817", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.12524751", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.12524751", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.53792817", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.53792817", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.12524751", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06737416", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.28936672", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.53792817", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.53792817", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.28936672", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06737416", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.01568694", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06737416", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.12524751", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.12524751", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06737416", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.01568694", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "0" 
     .YSet "6", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "1", "0", "0" 
     .YSet "5", "wavelength_centre/2", "70" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "-275.56759606" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "-275.56759606" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "-275.56759606" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "-275.56759606" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "-275.56759606" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "-275.56759606" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "-165.34055764" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "-165.34055764" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "-165.34055764" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "-165.34055764" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "-165.34055764" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "-165.34055764" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "-55.11351921" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "-55.11351921" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-55.11351921" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-55.11351921" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "-55.11351921" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "-55.11351921" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "55.11351921" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "55.11351921" 
     .Antenna "-11.53846", "11.53846", "0", "1", "55.11351921" 
     .Antenna "11.53846", "11.53846", "0", "1", "55.11351921" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "55.11351921" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "55.11351921" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "165.34055764" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "165.34055764" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "165.34055764" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "165.34055764" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "165.34055764" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "165.34055764" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "275.56759606" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "275.56759606" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "275.56759606" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "275.56759606" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "275.56759606" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "275.56759606" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "110" 
     .YSet "6", "wavelength_centre/2", "110" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle2" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "-220.45407685" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "-330.68111528" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "-440.9081537" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "-551.13519213" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "-661.36223055" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "-771.58926898" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "-330.68111528" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "-440.9081537" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "-551.13519213" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "-661.36223055" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "-771.58926898" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "-881.8163074" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "-440.9081537" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "-551.13519213" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-661.36223055" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-771.58926898" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "-881.8163074" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "-992.04334583" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "-551.13519213" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "-661.36223055" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-771.58926898" 
     .Antenna "11.53846", "11.53846", "0", "1", "-881.8163074" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "-992.04334583" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "-1102.27038425" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "-661.36223055" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "-771.58926898" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "-881.8163074" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "-992.04334583" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "-1102.27038425" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "-1212.49742268" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "-771.58926898" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "-881.8163074" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "-992.04334583" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "-1102.27038425" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "-1212.49742268" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "-1322.7244611" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "220.45407685" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "330.68111528" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "440.9081537" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "551.13519213" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "661.36223055" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "771.58926898" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "330.68111528" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "440.9081537" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "551.13519213" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "661.36223055" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "771.58926898" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "881.8163074" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "440.9081537" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "551.13519213" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "661.36223055" 
     .Antenna "11.53846", "-11.53846", "0", "1", "771.58926898" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "881.8163074" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "992.04334583" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "551.13519213" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "661.36223055" 
     .Antenna "-11.53846", "11.53846", "0", "1", "771.58926898" 
     .Antenna "11.53846", "11.53846", "0", "1", "881.8163074" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "992.04334583" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "1102.27038425" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "661.36223055" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "771.58926898" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "881.8163074" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "992.04334583" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "1102.27038425" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "1212.49742268" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "771.58926898" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "881.8163074" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "992.04334583" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "1102.27038425" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "1212.49742268" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "1322.7244611" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-220.45407685" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-330.68111528" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-440.9081537" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-551.13519213" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-330.68111528" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-440.9081537" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-551.13519213" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-661.36223055" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-440.9081537" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-551.13519213" 
     .Antenna "11.53846", "11.53846", "0", "1", "-661.36223055" 
     .Antenna "34.61538", "11.53846", "0", "1", "-771.58926898" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-551.13519213" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-661.36223055" 
     .Antenna "11.53846", "34.61538", "0", "1", "-771.58926898" 
     .Antenna "34.61538", "34.61538", "0", "1", "-881.8163074" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-0" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-0" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-0" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-0" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-0" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-0" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-0" 
     .Antenna "11.53846", "11.53846", "0", "1", "-0" 
     .Antenna "34.61538", "11.53846", "0", "1", "-0" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-0" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-0" 
     .Antenna "11.53846", "34.61538", "0", "1", "-0" 
     .Antenna "34.61538", "34.61538", "0", "1", "-0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-220.45407685" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-330.68111528" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-440.9081537" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-551.13519213" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-330.68111528" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-440.9081537" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-551.13519213" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-661.36223055" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-440.9081537" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-551.13519213" 
     .Antenna "11.53846", "11.53846", "0", "1", "-661.36223055" 
     .Antenna "34.61538", "11.53846", "0", "1", "-771.58926898" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-551.13519213" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-661.36223055" 
     .Antenna "11.53846", "34.61538", "0", "1", "-771.58926898" 
     .Antenna "34.61538", "34.61538", "0", "1", "-881.8163074" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-130.45407685" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-240.68111528" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-350.9081537" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-461.13519213" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-150.68111528" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-260.9081537" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-371.13519213" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-481.36223055" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-170.9081537" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-281.13519213" 
     .Antenna "11.53846", "11.53846", "0", "1", "-391.36223055" 
     .Antenna "34.61538", "11.53846", "0", "1", "-501.58926898" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-191.13519213" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-301.36223055" 
     .Antenna "11.53846", "34.61538", "0", "1", "-411.58926898" 
     .Antenna "34.61538", "34.61538", "0", "1", "-521.8163074" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-220.45407685" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-330.68111528" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-440.9081537" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-551.13519213" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-330.68111528" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-440.9081537" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-551.13519213" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-661.36223055" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-440.9081537" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-551.13519213" 
     .Antenna "11.53846", "11.53846", "0", "1", "-661.36223055" 
     .Antenna "34.61538", "11.53846", "0", "1", "-771.58926898" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-551.13519213" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-661.36223055" 
     .Antenna "11.53846", "34.61538", "0", "1", "-771.58926898" 
     .Antenna "34.61538", "34.61538", "0", "1", "-881.8163074" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.652395769231", "-57.652395769231", "0", "1", "191.135" 
     .Antenna "-34.591437461538", "-57.652395769231", "0", "1", "80.908" 
     .Antenna "-11.530479153846", "-57.652395769231", "0", "1", "330.681" 
     .Antenna "11.530479153846", "-57.652395769231", "0", "1", "220.454" 
     .Antenna "34.591437461538", "-57.652395769231", "0", "1", "110.227" 
     .Antenna "57.652395769231", "-57.652395769231", "0", "1", "0" 
     .Antenna "-57.652395769231", "-34.591437461538", "0", "1", "80.908" 
     .Antenna "-34.591437461538", "-34.591437461538", "0", "1", "330.681" 
     .Antenna "-11.530479153846", "-34.591437461538", "0", "1", "220.454" 
     .Antenna "11.530479153846", "-34.591437461538", "0", "1", "110.227" 
     .Antenna "34.591437461538", "-34.591437461538", "0", "1", "0" 
     .Antenna "57.652395769231", "-34.591437461538", "0", "1", "-110.227" 
     .Antenna "-57.652395769231", "-11.530479153846", "0", "1", "330.681" 
     .Antenna "-34.591437461538", "-11.530479153846", "0", "1", "220.454" 
     .Antenna "-11.530479153846", "-11.530479153846", "0", "1", "110.227" 
     .Antenna "11.530479153846", "-11.530479153846", "0", "1", "0" 
     .Antenna "34.591437461538", "-11.530479153846", "0", "1", "-110.227" 
     .Antenna "57.652395769231", "-11.530479153846", "0", "1", "-220.454" 
     .Antenna "-57.652395769231", "11.530479153846", "0", "1", "220.454" 
     .Antenna "-34.591437461538", "11.530479153846", "0", "1", "110.227" 
     .Antenna "-11.530479153846", "11.530479153846", "0", "1", "0" 
     .Antenna "11.530479153846", "11.530479153846", "0", "1", "-110.227" 
     .Antenna "34.591437461538", "11.530479153846", "0", "1", "-220.454" 
     .Antenna "57.652395769231", "11.530479153846", "0", "1", "-330.681" 
     .Antenna "-57.652395769231", "34.591437461538", "0", "1", "110.227" 
     .Antenna "-34.591437461538", "34.591437461538", "0", "1", "0" 
     .Antenna "-11.530479153846", "34.591437461538", "0", "1", "-110.227" 
     .Antenna "11.530479153846", "34.591437461538", "0", "1", "-220.454" 
     .Antenna "34.591437461538", "34.591437461538", "0", "1", "-330.681" 
     .Antenna "57.652395769231", "34.591437461538", "0", "1", "-80.908" 
     .Antenna "-57.652395769231", "57.652395769231", "0", "1", "0" 
     .Antenna "-34.591437461538", "57.652395769231", "0", "1", "-110.227" 
     .Antenna "-11.530479153846", "57.652395769231", "0", "1", "-220.454" 
     .Antenna "11.530479153846", "57.652395769231", "0", "1", "-330.681" 
     .Antenna "34.591437461538", "57.652395769231", "0", "1", "-80.908" 
     .Antenna "57.652395769231", "57.652395769231", "0", "1", "-191.135" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.652395769231", "-57.652395769231", "0", "1", "0" 
     .Antenna "-34.591437461538", "-57.652395769231", "0", "1", "-110.227" 
     .Antenna "-11.530479153846", "-57.652395769231", "0", "1", "-220.454" 
     .Antenna "11.530479153846", "-57.652395769231", "0", "1", "-330.681" 
     .Antenna "34.591437461538", "-57.652395769231", "0", "1", "-80.908" 
     .Antenna "57.652395769231", "-57.652395769231", "0", "1", "-191.135" 
     .Antenna "-57.652395769231", "-34.591437461538", "0", "1", "110.227" 
     .Antenna "-34.591437461538", "-34.591437461538", "0", "1", "0" 
     .Antenna "-11.530479153846", "-34.591437461538", "0", "1", "-110.227" 
     .Antenna "11.530479153846", "-34.591437461538", "0", "1", "-220.454" 
     .Antenna "34.591437461538", "-34.591437461538", "0", "1", "-330.681" 
     .Antenna "57.652395769231", "-34.591437461538", "0", "1", "-80.908" 
     .Antenna "-57.652395769231", "-11.530479153846", "0", "1", "220.454" 
     .Antenna "-34.591437461538", "-11.530479153846", "0", "1", "110.227" 
     .Antenna "-11.530479153846", "-11.530479153846", "0", "1", "0" 
     .Antenna "11.530479153846", "-11.530479153846", "0", "1", "-110.227" 
     .Antenna "34.591437461538", "-11.530479153846", "0", "1", "-220.454" 
     .Antenna "57.652395769231", "-11.530479153846", "0", "1", "-330.681" 
     .Antenna "-57.652395769231", "11.530479153846", "0", "1", "330.681" 
     .Antenna "-34.591437461538", "11.530479153846", "0", "1", "220.454" 
     .Antenna "-11.530479153846", "11.530479153846", "0", "1", "110.227" 
     .Antenna "11.530479153846", "11.530479153846", "0", "1", "0" 
     .Antenna "34.591437461538", "11.530479153846", "0", "1", "-110.227" 
     .Antenna "57.652395769231", "11.530479153846", "0", "1", "-220.454" 
     .Antenna "-57.652395769231", "34.591437461538", "0", "1", "80.908" 
     .Antenna "-34.591437461538", "34.591437461538", "0", "1", "330.681" 
     .Antenna "-11.530479153846", "34.591437461538", "0", "1", "220.454" 
     .Antenna "11.530479153846", "34.591437461538", "0", "1", "110.227" 
     .Antenna "34.591437461538", "34.591437461538", "0", "1", "0" 
     .Antenna "57.652395769231", "34.591437461538", "0", "1", "-110.227" 
     .Antenna "-57.652395769231", "57.652395769231", "0", "1", "191.135" 
     .Antenna "-34.591437461538", "57.652395769231", "0", "1", "80.908" 
     .Antenna "-11.530479153846", "57.652395769231", "0", "1", "330.681" 
     .Antenna "11.530479153846", "57.652395769231", "0", "1", "220.454" 
     .Antenna "34.591437461538", "57.652395769231", "0", "1", "110.227" 
     .Antenna "57.652395769231", "57.652395769231", "0", "1", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "139.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "29.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "279.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "1", "168.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "1", "58.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "1", "308.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "29.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "279.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "168.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "1", "58.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "1", "308.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "1", "198.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "279.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "168.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "58.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "308.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "1", "198.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "1", "87.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "1", "168.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "1", "58.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "308.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "198.1836926" 
     .Antenna "34.61538", "11.53846", "0", "1", "87.95665417" 
     .Antenna "57.69231", "11.53846", "0", "1", "337.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "1", "58.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "1", "308.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "1", "198.1836926" 
     .Antenna "11.53846", "34.61538", "0", "1", "87.95665417" 
     .Antenna "34.61538", "34.61538", "0", "1", "337.72961575" 
     .Antenna "57.69231", "34.61538", "0", "1", "227.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "1", "308.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "1", "198.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "1", "87.95665417" 
     .Antenna "11.53846", "57.69231", "0", "1", "337.72961575" 
     .Antenna "34.61538", "57.69231", "0", "1", "227.50257732" 
     .Antenna "57.69231", "57.69231", "0", "1", "117.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "-139.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "-29.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "-279.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "1", "-168.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "1", "-58.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "1", "-308.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "-29.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "-279.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "-168.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "1", "-58.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "1", "-308.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "1", "-198.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "-279.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "-168.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-58.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-308.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "1", "-198.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "1", "-87.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "1", "-168.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "1", "-58.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-308.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "-198.1836926" 
     .Antenna "34.61538", "11.53846", "0", "1", "-87.95665417" 
     .Antenna "57.69231", "11.53846", "0", "1", "-337.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "1", "-58.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "1", "-308.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "1", "-198.1836926" 
     .Antenna "11.53846", "34.61538", "0", "1", "-87.95665417" 
     .Antenna "34.61538", "34.61538", "0", "1", "-337.72961575" 
     .Antenna "57.69231", "34.61538", "0", "1", "-227.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "1", "-308.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "1", "-198.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "1", "-87.95665417" 
     .Antenna "11.53846", "57.69231", "0", "1", "-337.72961575" 
     .Antenna "34.61538", "57.69231", "0", "1", "-227.50257732" 
     .Antenna "57.69231", "57.69231", "0", "1", "-117.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "-110.2270" 
     .YSet "6", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "6", "wavelength_centre/2", "-110.34" 
     .YSet "6", "wavelength_centre/2", "-110.34" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-110.34" 
     .YSet "4", "wavelength_centre/2", "-110.34" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-63.703" 
     .YSet "4", "wavelength_centre/2", "-63.703" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "-220.45407685" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "-330.68111528" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "-440.9081537" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "-551.13519213" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "-661.36223055" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "-771.58926898" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "-330.68111528" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "-440.9081537" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "-551.13519213" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "-661.36223055" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "-771.58926898" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "-881.8163074" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "-440.9081537" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "-551.13519213" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "-661.36223055" 
     .Antenna "11.53846", "-11.53846", "0", "1", "-771.58926898" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "-881.8163074" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "-992.04334583" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "-551.13519213" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "-661.36223055" 
     .Antenna "-11.53846", "11.53846", "0", "1", "-771.58926898" 
     .Antenna "11.53846", "11.53846", "0", "1", "-881.8163074" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "-992.04334583" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "-1102.27038425" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "-661.36223055" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "-771.58926898" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "-881.8163074" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "-992.04334583" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "-1102.27038425" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "-1212.49742268" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "-771.58926898" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "-881.8163074" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "-992.04334583" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "-1102.27038425" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "-1212.49742268" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "-1322.7244611" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "139.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "29.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "279.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "1", "168.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "1", "58.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "1", "308.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "29.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "279.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "168.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "1", "58.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "1", "308.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "1", "198.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "279.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "168.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "58.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "308.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "1", "198.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "1", "87.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "1", "168.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "1", "58.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "308.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "198.1836926" 
     .Antenna "34.61538", "11.53846", "0", "1", "87.95665417" 
     .Antenna "57.69231", "11.53846", "0", "1", "337.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "1", "58.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "1", "308.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "1", "198.1836926" 
     .Antenna "11.53846", "34.61538", "0", "1", "87.95665417" 
     .Antenna "34.61538", "34.61538", "0", "1", "337.72961575" 
     .Antenna "57.69231", "34.61538", "0", "1", "227.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "1", "308.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "1", "198.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "1", "87.95665417" 
     .Antenna "11.53846", "57.69231", "0", "1", "337.72961575" 
     .Antenna "34.61538", "57.69231", "0", "1", "227.50257732" 
     .Antenna "57.69231", "57.69231", "0", "1", "117.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "252.34628196" 
     .Antenna "11.53846", "-57.69231", "0", "1", "96.46170928" 
     .Antenna "34.61538", "-57.69231", "0", "1", "300.57713659" 
     .Antenna "57.69231", "-57.69231", "0", "1", "144.69256391" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "252.34628196" 
     .Antenna "11.53846", "-34.61538", "0", "1", "96.46170928" 
     .Antenna "34.61538", "-34.61538", "0", "1", "300.57713659" 
     .Antenna "57.69231", "-34.61538", "0", "1", "144.69256391" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "252.34628196" 
     .Antenna "11.53846", "-11.53846", "0", "1", "96.46170928" 
     .Antenna "34.61538", "-11.53846", "0", "1", "300.57713659" 
     .Antenna "57.69231", "-11.53846", "0", "1", "144.69256391" 
     .Antenna "-57.69231", "11.53846", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "11.53846", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "11.53846", "0", "1", "252.34628196" 
     .Antenna "11.53846", "11.53846", "0", "1", "96.46170928" 
     .Antenna "34.61538", "11.53846", "0", "1", "300.57713659" 
     .Antenna "57.69231", "11.53846", "0", "1", "144.69256391" 
     .Antenna "-57.69231", "34.61538", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "34.61538", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "34.61538", "0", "1", "252.34628196" 
     .Antenna "11.53846", "34.61538", "0", "1", "96.46170928" 
     .Antenna "34.61538", "34.61538", "0", "1", "300.57713659" 
     .Antenna "57.69231", "34.61538", "0", "1", "144.69256391" 
     .Antenna "-57.69231", "57.69231", "0", "1", "204.11542732" 
     .Antenna "-34.61538", "57.69231", "0", "1", "48.23085464" 
     .Antenna "-11.53846", "57.69231", "0", "1", "252.34628196" 
     .Antenna "11.53846", "57.69231", "0", "1", "96.46170928" 
     .Antenna "34.61538", "57.69231", "0", "1", "300.57713659" 
     .Antenna "57.69231", "57.69231", "0", "1", "144.69256391" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-103.84615", "-103.84615", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "-103.84615", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "-103.84615", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "-103.84615", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "-103.84615", "0", "1", "300.57713659" 
     .Antenna "11.53846", "-103.84615", "0", "1", "144.69256391" 
     .Antenna "34.61538", "-103.84615", "0", "1", "348.80799123" 
     .Antenna "57.69231", "-103.84615", "0", "1", "192.92341855" 
     .Antenna "80.76923", "-103.84615", "0", "1", "37.03884587" 
     .Antenna "103.84615", "-103.84615", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "-80.76923", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "-80.76923", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "-80.76923", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "-80.76923", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "-80.76923", "0", "1", "300.57713659" 
     .Antenna "11.53846", "-80.76923", "0", "1", "144.69256391" 
     .Antenna "34.61538", "-80.76923", "0", "1", "348.80799123" 
     .Antenna "57.69231", "-80.76923", "0", "1", "192.92341855" 
     .Antenna "80.76923", "-80.76923", "0", "1", "37.03884587" 
     .Antenna "103.84615", "-80.76923", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "-57.69231", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "-57.69231", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "300.57713659" 
     .Antenna "11.53846", "-57.69231", "0", "1", "144.69256391" 
     .Antenna "34.61538", "-57.69231", "0", "1", "348.80799123" 
     .Antenna "57.69231", "-57.69231", "0", "1", "192.92341855" 
     .Antenna "80.76923", "-57.69231", "0", "1", "37.03884587" 
     .Antenna "103.84615", "-57.69231", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "-34.61538", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "-34.61538", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "300.57713659" 
     .Antenna "11.53846", "-34.61538", "0", "1", "144.69256391" 
     .Antenna "34.61538", "-34.61538", "0", "1", "348.80799123" 
     .Antenna "57.69231", "-34.61538", "0", "1", "192.92341855" 
     .Antenna "80.76923", "-34.61538", "0", "1", "37.03884587" 
     .Antenna "103.84615", "-34.61538", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "-11.53846", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "-11.53846", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "300.57713659" 
     .Antenna "11.53846", "-11.53846", "0", "1", "144.69256391" 
     .Antenna "34.61538", "-11.53846", "0", "1", "348.80799123" 
     .Antenna "57.69231", "-11.53846", "0", "1", "192.92341855" 
     .Antenna "80.76923", "-11.53846", "0", "1", "37.03884587" 
     .Antenna "103.84615", "-11.53846", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "11.53846", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "11.53846", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "11.53846", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "11.53846", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "11.53846", "0", "1", "300.57713659" 
     .Antenna "11.53846", "11.53846", "0", "1", "144.69256391" 
     .Antenna "34.61538", "11.53846", "0", "1", "348.80799123" 
     .Antenna "57.69231", "11.53846", "0", "1", "192.92341855" 
     .Antenna "80.76923", "11.53846", "0", "1", "37.03884587" 
     .Antenna "103.84615", "11.53846", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "34.61538", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "34.61538", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "34.61538", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "34.61538", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "34.61538", "0", "1", "300.57713659" 
     .Antenna "11.53846", "34.61538", "0", "1", "144.69256391" 
     .Antenna "34.61538", "34.61538", "0", "1", "348.80799123" 
     .Antenna "57.69231", "34.61538", "0", "1", "192.92341855" 
     .Antenna "80.76923", "34.61538", "0", "1", "37.03884587" 
     .Antenna "103.84615", "34.61538", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "57.69231", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "57.69231", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "57.69231", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "57.69231", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "57.69231", "0", "1", "300.57713659" 
     .Antenna "11.53846", "57.69231", "0", "1", "144.69256391" 
     .Antenna "34.61538", "57.69231", "0", "1", "348.80799123" 
     .Antenna "57.69231", "57.69231", "0", "1", "192.92341855" 
     .Antenna "80.76923", "57.69231", "0", "1", "37.03884587" 
     .Antenna "103.84615", "57.69231", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "80.76923", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "80.76923", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "80.76923", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "80.76923", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "80.76923", "0", "1", "300.57713659" 
     .Antenna "11.53846", "80.76923", "0", "1", "144.69256391" 
     .Antenna "34.61538", "80.76923", "0", "1", "348.80799123" 
     .Antenna "57.69231", "80.76923", "0", "1", "192.92341855" 
     .Antenna "80.76923", "80.76923", "0", "1", "37.03884587" 
     .Antenna "103.84615", "80.76923", "0", "1", "241.15427319" 
     .Antenna "-103.84615", "103.84615", "0", "1", "204.11542732" 
     .Antenna "-80.76923", "103.84615", "0", "1", "48.23085464" 
     .Antenna "-57.69231", "103.84615", "0", "1", "252.34628196" 
     .Antenna "-34.61538", "103.84615", "0", "1", "96.46170928" 
     .Antenna "-11.53846", "103.84615", "0", "1", "300.57713659" 
     .Antenna "11.53846", "103.84615", "0", "1", "144.69256391" 
     .Antenna "34.61538", "103.84615", "0", "1", "348.80799123" 
     .Antenna "57.69231", "103.84615", "0", "1", "192.92341855" 
     .Antenna "80.76923", "103.84615", "0", "1", "37.03884587" 
     .Antenna "103.84615", "103.84615", "0", "1", "241.15427319" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle2" 
     .Theta "60" 
     .Phi "60" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-103.84615", "-103.84615", "0", "1", "139.54592315" 
     .Antenna "-80.76923", "-103.84615", "0", "1", "29.31888472" 
     .Antenna "-57.69231", "-103.84615", "0", "1", "279.0918463" 
     .Antenna "-34.61538", "-103.84615", "0", "1", "168.86480787" 
     .Antenna "-11.53846", "-103.84615", "0", "1", "58.63776945" 
     .Antenna "11.53846", "-103.84615", "0", "1", "308.41073102" 
     .Antenna "34.61538", "-103.84615", "0", "1", "198.1836926" 
     .Antenna "57.69231", "-103.84615", "0", "1", "87.95665417" 
     .Antenna "80.76923", "-103.84615", "0", "1", "337.72961575" 
     .Antenna "103.84615", "-103.84615", "0", "1", "227.50257732" 
     .Antenna "-103.84615", "-80.76923", "0", "1", "29.31888472" 
     .Antenna "-80.76923", "-80.76923", "0", "1", "279.0918463" 
     .Antenna "-57.69231", "-80.76923", "0", "1", "168.86480787" 
     .Antenna "-34.61538", "-80.76923", "0", "1", "58.63776945" 
     .Antenna "-11.53846", "-80.76923", "0", "1", "308.41073102" 
     .Antenna "11.53846", "-80.76923", "0", "1", "198.1836926" 
     .Antenna "34.61538", "-80.76923", "0", "1", "87.95665417" 
     .Antenna "57.69231", "-80.76923", "0", "1", "337.72961575" 
     .Antenna "80.76923", "-80.76923", "0", "1", "227.50257732" 
     .Antenna "103.84615", "-80.76923", "0", "1", "117.2755389" 
     .Antenna "-103.84615", "-57.69231", "0", "1", "279.0918463" 
     .Antenna "-80.76923", "-57.69231", "0", "1", "168.86480787" 
     .Antenna "-57.69231", "-57.69231", "0", "1", "58.63776945" 
     .Antenna "-34.61538", "-57.69231", "0", "1", "308.41073102" 
     .Antenna "-11.53846", "-57.69231", "0", "1", "198.1836926" 
     .Antenna "11.53846", "-57.69231", "0", "1", "87.95665417" 
     .Antenna "34.61538", "-57.69231", "0", "1", "337.72961575" 
     .Antenna "57.69231", "-57.69231", "0", "1", "227.50257732" 
     .Antenna "80.76923", "-57.69231", "0", "1", "117.2755389" 
     .Antenna "103.84615", "-57.69231", "0", "1", "7.04850047" 
     .Antenna "-103.84615", "-34.61538", "0", "1", "168.86480787" 
     .Antenna "-80.76923", "-34.61538", "0", "1", "58.63776945" 
     .Antenna "-57.69231", "-34.61538", "0", "1", "308.41073102" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "198.1836926" 
     .Antenna "-11.53846", "-34.61538", "0", "1", "87.95665417" 
     .Antenna "11.53846", "-34.61538", "0", "1", "337.72961575" 
     .Antenna "34.61538", "-34.61538", "0", "1", "227.50257732" 
     .Antenna "57.69231", "-34.61538", "0", "1", "117.2755389" 
     .Antenna "80.76923", "-34.61538", "0", "1", "7.04850047" 
     .Antenna "103.84615", "-34.61538", "0", "1", "256.82146205" 
     .Antenna "-103.84615", "-11.53846", "0", "1", "58.63776945" 
     .Antenna "-80.76923", "-11.53846", "0", "1", "308.41073102" 
     .Antenna "-57.69231", "-11.53846", "0", "1", "198.1836926" 
     .Antenna "-34.61538", "-11.53846", "0", "1", "87.95665417" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "337.72961575" 
     .Antenna "11.53846", "-11.53846", "0", "1", "227.50257732" 
     .Antenna "34.61538", "-11.53846", "0", "1", "117.2755389" 
     .Antenna "57.69231", "-11.53846", "0", "1", "7.04850047" 
     .Antenna "80.76923", "-11.53846", "0", "1", "256.82146205" 
     .Antenna "103.84615", "-11.53846", "0", "1", "146.59442362" 
     .Antenna "-103.84615", "11.53846", "0", "1", "308.41073102" 
     .Antenna "-80.76923", "11.53846", "0", "1", "198.1836926" 
     .Antenna "-57.69231", "11.53846", "0", "1", "87.95665417" 
     .Antenna "-34.61538", "11.53846", "0", "1", "337.72961575" 
     .Antenna "-11.53846", "11.53846", "0", "1", "227.50257732" 
     .Antenna "11.53846", "11.53846", "0", "1", "117.2755389" 
     .Antenna "34.61538", "11.53846", "0", "1", "7.04850047" 
     .Antenna "57.69231", "11.53846", "0", "1", "256.82146205" 
     .Antenna "80.76923", "11.53846", "0", "1", "146.59442362" 
     .Antenna "103.84615", "11.53846", "0", "1", "36.3673852" 
     .Antenna "-103.84615", "34.61538", "0", "1", "198.1836926" 
     .Antenna "-80.76923", "34.61538", "0", "1", "87.95665417" 
     .Antenna "-57.69231", "34.61538", "0", "1", "337.72961575" 
     .Antenna "-34.61538", "34.61538", "0", "1", "227.50257732" 
     .Antenna "-11.53846", "34.61538", "0", "1", "117.2755389" 
     .Antenna "11.53846", "34.61538", "0", "1", "7.04850047" 
     .Antenna "34.61538", "34.61538", "0", "1", "256.82146205" 
     .Antenna "57.69231", "34.61538", "0", "1", "146.59442362" 
     .Antenna "80.76923", "34.61538", "0", "1", "36.3673852" 
     .Antenna "103.84615", "34.61538", "0", "1", "286.14034677" 
     .Antenna "-103.84615", "57.69231", "0", "1", "87.95665417" 
     .Antenna "-80.76923", "57.69231", "0", "1", "337.72961575" 
     .Antenna "-57.69231", "57.69231", "0", "1", "227.50257732" 
     .Antenna "-34.61538", "57.69231", "0", "1", "117.2755389" 
     .Antenna "-11.53846", "57.69231", "0", "1", "7.04850047" 
     .Antenna "11.53846", "57.69231", "0", "1", "256.82146205" 
     .Antenna "34.61538", "57.69231", "0", "1", "146.59442362" 
     .Antenna "57.69231", "57.69231", "0", "1", "36.3673852" 
     .Antenna "80.76923", "57.69231", "0", "1", "286.14034677" 
     .Antenna "103.84615", "57.69231", "0", "1", "175.91330835" 
     .Antenna "-103.84615", "80.76923", "0", "1", "337.72961575" 
     .Antenna "-80.76923", "80.76923", "0", "1", "227.50257732" 
     .Antenna "-57.69231", "80.76923", "0", "1", "117.2755389" 
     .Antenna "-34.61538", "80.76923", "0", "1", "7.04850047" 
     .Antenna "-11.53846", "80.76923", "0", "1", "256.82146205" 
     .Antenna "11.53846", "80.76923", "0", "1", "146.59442362" 
     .Antenna "34.61538", "80.76923", "0", "1", "36.3673852" 
     .Antenna "57.69231", "80.76923", "0", "1", "286.14034677" 
     .Antenna "80.76923", "80.76923", "0", "1", "175.91330835" 
     .Antenna "103.84615", "80.76923", "0", "1", "65.68626992" 
     .Antenna "-103.84615", "103.84615", "0", "1", "227.50257732" 
     .Antenna "-80.76923", "103.84615", "0", "1", "117.2755389" 
     .Antenna "-57.69231", "103.84615", "0", "1", "7.04850047" 
     .Antenna "-34.61538", "103.84615", "0", "1", "256.82146205" 
     .Antenna "-11.53846", "103.84615", "0", "1", "146.59442362" 
     .Antenna "11.53846", "103.84615", "0", "1", "36.3673852" 
     .Antenna "34.61538", "103.84615", "0", "1", "286.14034677" 
     .Antenna "57.69231", "103.84615", "0", "1", "175.91330835" 
     .Antenna "80.76923", "103.84615", "0", "1", "65.68626992" 
     .Antenna "103.84615", "103.84615", "0", "1", "315.4592315" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "139.54592315" 
     .Antenna "-11.53846", "-34.61538", "0", "0.41903175", "29.31888472" 
     .Antenna "11.53846", "-34.61538", "0", "0.41903175", "279.0918463" 
     .Antenna "34.61538", "-34.61538", "0", "1", "168.86480787" 
     .Antenna "-34.61538", "-11.53846", "0", "0.41903175", "29.31888472" 
     .Antenna "-11.53846", "-11.53846", "0", "0.17558761", "279.0918463" 
     .Antenna "11.53846", "-11.53846", "0", "0.17558761", "168.86480787" 
     .Antenna "34.61538", "-11.53846", "0", "0.41903175", "58.63776945" 
     .Antenna "-34.61538", "11.53846", "0", "0.41903175", "279.0918463" 
     .Antenna "-11.53846", "11.53846", "0", "0.17558761", "168.86480787" 
     .Antenna "11.53846", "11.53846", "0", "0.17558761", "58.63776945" 
     .Antenna "34.61538", "11.53846", "0", "0.41903175", "308.41073102" 
     .Antenna "-34.61538", "34.61538", "0", "1", "168.86480787" 
     .Antenna "-11.53846", "34.61538", "0", "0.41903175", "58.63776945" 
     .Antenna "11.53846", "34.61538", "0", "0.41903175", "308.41073102" 
     .Antenna "34.61538", "34.61538", "0", "1", "198.1836926" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "FALSE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "1", "180" 
     .Antenna "-11.53846", "-34.61538", "0", "0.41903175", "360" 
     .Antenna "11.53846", "-34.61538", "0", "0.41903175", "180" 
     .Antenna "34.61538", "-34.61538", "0", "1", "360" 
     .Antenna "-34.61538", "-11.53846", "0", "0.41903175", "180" 
     .Antenna "-11.53846", "-11.53846", "0", "0.17558761", "360" 
     .Antenna "11.53846", "-11.53846", "0", "0.17558761", "180" 
     .Antenna "34.61538", "-11.53846", "0", "0.41903175", "360" 
     .Antenna "-34.61538", "11.53846", "0", "0.41903175", "180" 
     .Antenna "-11.53846", "11.53846", "0", "0.17558761", "360" 
     .Antenna "11.53846", "11.53846", "0", "0.17558761", "180" 
     .Antenna "34.61538", "11.53846", "0", "0.41903175", "360" 
     .Antenna "-34.61538", "34.61538", "0", "1", "180" 
     .Antenna "-11.53846", "34.61538", "0", "0.41903175", "360" 
     .Antenna "11.53846", "34.61538", "0", "0.41903175", "180" 
     .Antenna "34.61538", "34.61538", "0", "1", "360" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "0" 
     .YSet "4", "wavelength_centre/2", "0" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-120.2270" 
     .YSet "4", "wavelength_centre/2", "-120.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle2" 
     .Theta "60" 
     .Phi "60" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-125.2270" 
     .YSet "4", "wavelength_centre/2", "-125.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "60" 
     .Phi "60" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "60" 
     .Phi "60" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-110.2270" 
     .YSet "4", "wavelength_centre/2", "-110.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-170.2270" 
     .YSet "4", "wavelength_centre/2", "-170.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Directivity" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "0" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "11.53846", "-11.53846", "0", "1", "0" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "-11.53846", "11.53846", "0", "1", "0" 
     .Antenna "11.53846", "11.53846", "0", "1", "0" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "0" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "0" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "0" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "0" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "0" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "0" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "0" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "0" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "0" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Rectangular" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .SetList
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01431209", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06342272", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11963314", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11963314", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06342272", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01431209", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06342272", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28105207", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53014344", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.53014344", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.28105207", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06342272", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06342272", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.28105207", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.53014344", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.53014344", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.28105207", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06342272", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01431209", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06342272", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11963314", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11963314", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06342272", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01431209", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01419813", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06291772", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11868056", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11868056", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06291772", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01419813", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.0632607", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28033406", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52878907", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.52878907", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.28033406", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.0632607", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.0632607", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.28033406", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.52878907", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.52878907", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.28033406", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.0632607", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01419813", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06291772", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11868056", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11868056", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06291772", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01419813", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01379791", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06114419", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06114419", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01379791", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.062682", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27776963", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.27776963", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.062682", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.062682", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.27776963", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.27776963", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.062682", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01379791", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06114419", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06114419", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01379791", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01340132", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.05938673", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11202011", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11202011", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.05938673", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01340132", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06209314", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27516016", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.51902964", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.51902964", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.27516016", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06209314", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06209314", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.27516016", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.51902964", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.51902964", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.27516016", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06209314", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01340132", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.05938673", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11202011", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11202011", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.05938673", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01340132", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01379791", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.062682", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11963314", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11963314", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.062682", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01379791", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06114419", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27776963", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53014344", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.53014344", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.27776963", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06114419", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06114419", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.27776963", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.53014344", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.53014344", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.27776963", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06114419", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01379791", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.062682", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11963314", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11963314", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.062682", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01379791", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "149.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "44.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "299.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "193.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "88.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "343.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "44.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "299.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "193.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "88.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "343.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "238.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "299.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "193.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "88.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "343.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "238.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "132.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "193.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "88.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "343.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "238.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "132.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "27.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "88.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "343.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "238.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "132.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "27.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "282.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "343.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "238.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "132.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "27.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "282.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "177.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01431209", "139.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06342272", "29.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11963314", "279.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11963314", "168.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06342272", "58.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01431209", "308.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06342272", "29.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28105207", "279.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53014344", "168.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.53014344", "58.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.28105207", "308.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06342272", "198.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "279.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "168.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "58.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "308.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "198.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "87.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "168.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "58.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "308.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "198.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "87.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "337.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06342272", "58.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.28105207", "308.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.53014344", "198.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.53014344", "87.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.28105207", "337.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06342272", "227.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01431209", "308.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06342272", "198.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11963314", "87.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11963314", "337.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06342272", "227.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01431209", "117.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01431209", "159.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06342272", "59.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11963314", "319.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11963314", "218.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06342272", "118.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01431209", "18.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06342272", "59.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28105207", "319.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53014344", "218.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.53014344", "118.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.28105207", "18.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06342272", "278.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "319.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "218.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "118.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "18.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "278.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "177.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "218.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "118.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "18.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "278.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "177.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "77.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06342272", "118.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.28105207", "18.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.53014344", "278.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.53014344", "177.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.28105207", "77.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06342272", "337.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01431209", "18.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06342272", "278.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11963314", "177.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11963314", "77.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06342272", "337.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01431209", "237.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01431209", "119.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06342272", "359.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11963314", "239.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11963314", "118.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06342272", "358.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01431209", "238.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06342272", "359.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.28105207", "239.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.53014344", "118.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.53014344", "358.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.28105207", "238.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06342272", "118.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11963314", "239.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.53014344", "118.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "358.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "238.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.53014344", "118.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11963314", "357.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11963314", "118.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.53014344", "358.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "238.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "118.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.53014344", "357.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11963314", "237.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06342272", "358.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.28105207", "238.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.53014344", "118.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.53014344", "357.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.28105207", "237.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06342272", "117.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01431209", "238.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06342272", "118.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11963314", "357.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11963314", "237.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06342272", "117.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01431209", "357.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "90" 
     .Phi "90" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.0133022", "119.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.06043008", "359.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11533518", "239.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11533518", "118.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.06043008", "358.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.0133022", "238.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.06043008", "359.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.27452552", "239.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.52395183", "118.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.52395183", "358.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.27452552", "238.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.06043008", "118.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11533518", "239.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.52395183", "118.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "358.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "238.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.52395183", "118.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11533518", "357.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11533518", "118.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.52395183", "358.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "238.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "118.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.52395183", "357.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11533518", "237.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.06043008", "358.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.27452552", "238.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.52395183", "118.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.52395183", "357.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.27452552", "237.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.06043008", "117.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.0133022", "238.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.06043008", "118.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11533518", "357.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11533518", "237.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.06043008", "117.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.0133022", "357.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Pfield" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


'@ farfield array properties

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldArray
     .Reset 
     .UseArray "TRUE" 
     .Arraytype "Edit" 
     .XSet "4", "wavelength_centre/2", "-130.2270" 
     .YSet "4", "wavelength_centre/2", "-130.2270" 
     .ZSet "1", "0", "0" 
     .Antenna "-57.69231", "-57.69231", "0", "0.01254851", "109.54592315" 
     .Antenna "-34.61538", "-57.69231", "0", "0.05814176", "344.31888472" 
     .Antenna "-11.53846", "-57.69231", "0", "0.11202011", "219.0918463" 
     .Antenna "11.53846", "-57.69231", "0", "0.11202011", "93.86480787" 
     .Antenna "34.61538", "-57.69231", "0", "0.05814176", "328.63776945" 
     .Antenna "57.69231", "-57.69231", "0", "0.01254851", "203.41073102" 
     .Antenna "-57.69231", "-34.61538", "0", "0.05814176", "344.31888472" 
     .Antenna "-34.61538", "-34.61538", "0", "0.26939177", "219.0918463" 
     .Antenna "-11.53846", "-34.61538", "0", "0.51902964", "93.86480787" 
     .Antenna "11.53846", "-34.61538", "0", "0.51902964", "328.63776945" 
     .Antenna "34.61538", "-34.61538", "0", "0.26939177", "203.41073102" 
     .Antenna "57.69231", "-34.61538", "0", "0.05814176", "78.1836926" 
     .Antenna "-57.69231", "-11.53846", "0", "0.11202011", "219.0918463" 
     .Antenna "-34.61538", "-11.53846", "0", "0.51902964", "93.86480787" 
     .Antenna "-11.53846", "-11.53846", "0", "1", "328.63776945" 
     .Antenna "11.53846", "-11.53846", "0", "1", "203.41073102" 
     .Antenna "34.61538", "-11.53846", "0", "0.51902964", "78.1836926" 
     .Antenna "57.69231", "-11.53846", "0", "0.11202011", "312.95665417" 
     .Antenna "-57.69231", "11.53846", "0", "0.11202011", "93.86480787" 
     .Antenna "-34.61538", "11.53846", "0", "0.51902964", "328.63776945" 
     .Antenna "-11.53846", "11.53846", "0", "1", "203.41073102" 
     .Antenna "11.53846", "11.53846", "0", "1", "78.1836926" 
     .Antenna "34.61538", "11.53846", "0", "0.51902964", "312.95665417" 
     .Antenna "57.69231", "11.53846", "0", "0.11202011", "187.72961575" 
     .Antenna "-57.69231", "34.61538", "0", "0.05814176", "328.63776945" 
     .Antenna "-34.61538", "34.61538", "0", "0.26939177", "203.41073102" 
     .Antenna "-11.53846", "34.61538", "0", "0.51902964", "78.1836926" 
     .Antenna "11.53846", "34.61538", "0", "0.51902964", "312.95665417" 
     .Antenna "34.61538", "34.61538", "0", "0.26939177", "187.72961575" 
     .Antenna "57.69231", "34.61538", "0", "0.05814176", "62.50257732" 
     .Antenna "-57.69231", "57.69231", "0", "0.01254851", "203.41073102" 
     .Antenna "-34.61538", "57.69231", "0", "0.05814176", "78.1836926" 
     .Antenna "-11.53846", "57.69231", "0", "0.11202011", "312.95665417" 
     .Antenna "11.53846", "57.69231", "0", "0.11202011", "187.72961575" 
     .Antenna "34.61538", "57.69231", "0", "0.05814176", "62.50257732" 
     .Antenna "57.69231", "57.69231", "0", "0.01254851", "297.2755389" 
End With


'@ farfield plot options

'[VERSION]2018.0|27.0.2|20171026[/VERSION]
With FarfieldPlot 
     .Plottype "Polar" 
     .Vary "angle1" 
     .Theta "45" 
     .Phi "45" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "False" 
     .SymmetricRange "False" 
     .SetTimeDomainFF "False" 
     .SetFrequency "6.5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With 


