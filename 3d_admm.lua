PluginRequired("FluidOptim")
PluginRequired("PLaplacian")
PluginRequired("ADMMOptim")
--All the previous scripts are necessary without exception
-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("util/conv_rates_static.lua")
ug_load_script("util/solver_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("obstacle_optim_3d_util.lua")

ProfileLUA(true)
PrintBuildConfiguration()
--hlrn path
package.path =  package.path ..";/home/hhmpinzo/ug4/apps/plaplace_optim/lua-matrix/?.lua"
--hawk path
--package.path = package.path .. ";/zhome/academic/HLRS/xbo/xbopinzo/ug4/apps/plaplace_optim/lua-matrix/?.lua"
--package.path = package.path .. ";../../apps/plaplace_optim/lua-matrix/?.lua"
--package.path = package.path .. ";/afs/math.uni-hamburg.de/users/oa/bay8214/ug4/apps/plaplace_optim/lua-matrix/?.lua"
local matrix = require "matrix"
--[[
Try using the command:
mpirun -np 4 ugshell -ex ../../apps/plaplace_optim/ns_test.lua -visc 0.1 -eps 1.0e-8 -numRefs 1 -nsMaxIts 10 
-step_length 1.0 -control 1.0e-1 -nsTol 1.0e-10 -numSteps 1 -p_max 2.0
]]
--Constants
dim = 3
numPreRefs=0
vorder=2
porder=1
diameter=6
m=dim+1
flowNames="v1,v2,v3"
pressureName="p"
adjointFlowNames="q1,q2,q3"
adjointPressure="h"
flowOutputFile="vtk_flows"
adjointFlowOutputFile="vtk_adjointFlows"
nodalFlows="flows"
adjointNodalFlows="adjointFlows"


--Simulation Parameters
numRefs = util.GetParamNumber("-numRefs", 2)--refinements after grid partitioning/distribution
numSteps = util.GetParamNumber("-numSteps", 400)
admmSteps = util.GetParamNumber("-admmSteps", 2)
visc = util.GetParamNumber("-visc", 0.02)--medium viscosity
stab = util.GetParamNumber("-stab", 0.0)--stabilization for navier stokes
stabType = util.GetParamNumber("-stabType",0.0)
sigma_threshold = util.GetParamNumber("-sigma_threshold",0.3)--Threshold of gradient Nabla(u)
scaling = util.GetParamNumber("-scaling",1.0)--scaling for the sensitivity J'
admm_tolerance = util.GetParamNumber("-admm_tolerance", 1e-2)--tolerance for admm convergence
step_length=util.GetParamNumber("-step_length",1.0)--small number for Lu sensitivity
p_increase = util.GetParamNumber("-p_inc", 0.19)--increase interval for p
step_control = util.GetParamNumber("-control", 1.0)--control for p term
tau=util.GetParamNumber("-tau", 1e0)
norm_limit=util.GetParamNumber("-norm_limit", 1e-5)
bReadInitialGuess = util.GetParamNumber("-restart", -1)
gridName = util.GetParam("-grid", "./grids/box_3D_elongated.ugx")
--gridName = util.GetParam("-grid", "./grids/box_2D.ugx")
--Newton solver settings
nsMaxIts = util.GetParamNumber("-nsMaxIts", 30)
nsTol = util.GetParamNumber("-nsTol", 1e-9)

lambda_vol = util.GetParamNumber("-lambda_vol", 0.0)
lambda_x   = util.GetParamNumber("-lambda_x", 0.0)
lambda_y   = util.GetParamNumber("-lambda_y", 0.0)
lambda_z   = util.GetParamNumber("-lambda_z", 0.0)
--Boolean parameters
bNewtonOutput = util.GetParamBool("-bNewtonOutput", false)--output newton method per step convergence information
bOutputMesh  = util.GetParamBool("-bOutputMesh",true)--output VTK visualization files p-laplacian
bOutputFlows  = util.GetParamBool("-bOutputFlows",true)--output VTK visualization files flow and pressure
bOutputPressure  = util.GetParamBool("-bOutputPressure",false)--output VTK visualization files flow and pressure
bOutputAdjoints  = util.GetParamBool("-bOutputAdjoints",false)--output VTK visualization files of flow and pressureadjoints
bDebugOutput = util.GetParamBool("-bDebugOutput",false)--output VTK visualization files of Lu and RHS of big problem
bDebugNodalPositions = util.GetParamBool("-bDebugNodalPositions",false)
bDebugSensitivity = util.GetParamBool("-bDebugSensitivity",false)--output VTK files for J'
bDoNothing = util.GetParamBool("-bDoNothing",true)--do nothing on flow outlet
bOutputIntermediateUp=util.GetParamBool("-bOutputIntermediateUp",false)--output VTK files for intermediate u files
bActivateProfiler = util.GetParamBool("-bActivateProfiler",true)

print("THE PARAMETERS USED FOR EXECUTION ARE: ")
print("grid: "..gridName)
print("numPreRefs:   ".. numPreRefs)
print("numRefs:      ".. numRefs)
print("numSteps:     ".. numSteps)
print("press.grad. stabilization:".. stab)
print("average.based stab:".. stabType)
print("viscosity:".. visc)
print("velocity order "..vorder)
print("pressure order "..porder)
print("step_control:     ".. step_control)
print("step_length:     ".. step_length)
print("Newton solver tolerance set to"..nsTol)
print("restarted at step"..bReadInitialGuess)


step_control_init = step_control
step_length_init = step_length
p_current=2.0

-- initialize ug with the world dimension and the algebra type
InitUG(dim, AlgebraType("CPU", 1));

-- load grid into domain
dom = Domain()
LoadDomain(dom, gridName)
--dom = util.CreateDomain(gridName, 0, {})
print("Loaded domain from " .. gridName)
number_elements=dom:domain_info():num_surface_elements()
--[[
local refiner =  GlobalDomainRefiner(dom)
for i=1,numPreRefs do refiner:refine(); end
if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, false) == false then
	--print("Error while Distributing Grid. Aborting."); exit();
end
for i=numPreRefs+1,numRefs do refiner:refine(); end
--]]

----[[
-- This balancing setup makes sense for structured grids with uniform refinement
balancerDesc = {
	--[[
		partitioner = {
			name = "dynamicBisection",
			verbose = false,
			enableXCuts = true,
			enableYCuts = true,
			enableZCuts = true,
			longestSplitAxis = false,
			clusteredSiblings = true,
			balanceThreshold = 0.9,
			numSplitImprovements = 10
		},
	--]]
----[[
	partitioner =
	{
		name = "parmetis",
		balanceWeights = nil,
		communicationWeights = nil,
		options = nil,
		itrFactor               = 1000,
		verbose = false,
		clusteredSiblings = true,
		balanceThreshold = 0.9
	},
--]]
	hierarchy = {
		name                            = "standard",   --- ["standard", "lessRedists", "noRedists"]
		minElemsPerProcPerLevel         = 4,
		maxRedistProcs                  = 384,
		qualityRedistLevelOffset        = 1,
		intermediateRedistributions     = true,
		allowForRedistOnEachLevel       = false, -- DANGEROUS: 'true' only useful for non-adaptive runs with static partitioning. Leads to errors else.
		{-- level 0
			upperLvl = 0,
			maxProcs = 1
		},
		{-- levels 1
			upperLvl = 1,
			maxProcs = 96
		},
		{-- levels 2
			upperLvl = 2,
			maxProcs = 384
		},
		{-- levels 3
			upperLvl = 3,
			maxProcs = 1536
		},
		{-- levels 4
			upperLvl = 4,
			maxProcs = 6144
		},
		{-- levels 5
			upperLvl = 5,
			maxProcs = 12288
		}

	},
}

util.refinement.CreateRegularHierarchy(dom, numRefs, false, balancerDesc)
--]]

print(dom:domain_info():to_string())

--NAVIER STOKES 
--DIRICHLET BOUNDARY VALUES AT INLET
function InletVelocities(x,y,z,t)
	local s=math.sqrt(y*y+z*z)*math.pi/diameter
	return math.max(0.0, math.cos(s))
	--return 0
end
print("NAVIER STOKES: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
NavierStokes_ApproxSpace = ApproximationSpace(dom)
NavierStokes_ApproxSpace:add_fct(flowNames,"Lagrange", vorder)
NavierStokes_ApproxSpace:add_fct(pressureName,"Lagrange",porder) 
NavierStokes_ApproxSpace:init_levels()
NavierStokes_ApproxSpace:init_top_surface()

print("NAVIER STOKES: Approx. Space:")
NavierStokes_ApproxSpace:print_statistic()

NavierStokes_ElemDisc = IncompressibleNavierStokes("v1,v2,v3,p ", "outer")
--useful settings
NavierStokes_ElemDisc:set_nonlinear(true)
NavierStokes_ElemDisc:set_picard(false)
NavierStokes_ElemDisc:set_kinematic_viscosity(visc)
NavierStokes_ElemDisc:set_stabilization(stab)
NavierStokes_ElemDisc:set_stabilization_type(stabType)
--************BOUNDARY CONDITIONS*********--
NavierStokes_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
NavierStokes_Dirich:add("InletVelocities","v1","inlet")
NavierStokes_Dirich:add(0,"v2","inlet")
NavierStokes_Dirich:add(0,"v3","inlet")
--************WALL BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","wall")
NavierStokes_Dirich:add(0,"v2","wall")
NavierStokes_Dirich:add(0,"v3","wall")

--************OBSTACLE SURFACE BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","obstacle_surface")
NavierStokes_Dirich:add(0,"v2","obstacle_surface")
NavierStokes_Dirich:add(0,"v3","obstacle_surface")

if not bDoNothing then
	print("Output flows set")
	NavierStokes_Dirich:add("InletVelocities","v1","outlet")
	NavierStokes_Dirich:add(0,"v2","outlet")
    NavierStokes_Dirich:add(0,"v3","outlet")
end
-- Domain Discretization 
NavierStokes_DomainDisc = DomainDiscretization(NavierStokes_ApproxSpace)
NavierStokes_DomainDisc:add(NavierStokes_ElemDisc)
NavierStokes_DomainDisc:add(NavierStokes_Dirich)
print("NAVIER STOKES: create GridFunctions, Matrix Operator, and GlobalGridFunctionGradientDatas")
v = GridFunction(NavierStokes_ApproxSpace);v:set(0.0)
v_temp = GridFunction(NavierStokes_ApproxSpace);v_temp:set(0.0)
NavierStokes_DomainDisc:adjust_solution(v)
NavierStokes_DomainDisc:adjust_solution(v_temp)
v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
v1_value_global=GlobalGridFunctionNumberData(v,"v1")
v2_value_global=GlobalGridFunctionNumberData(v,"v2")
v3_value_global=GlobalGridFunctionNumberData(v,"v3")

---Pressure
p_value_global=GlobalGridFunctionNumberData(v,"p")
--Nonlinear Matrix Operator
navier_Op = AssembledOperator()
navier_Op:set_discretization(NavierStokes_DomainDisc)
print("NAVIER STOKES: all set")


print("ADJOINT FLOW: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
AdjointFlow_ApproxSpace = ApproximationSpace(dom)
AdjointFlow_ApproxSpace:add_fct(adjointFlowNames,"Lagrange",vorder)
AdjointFlow_ApproxSpace:add_fct(adjointPressure,"Lagrange",porder) 
AdjointFlow_ApproxSpace:init_levels()
AdjointFlow_ApproxSpace:init_top_surface()
print("ADJOINT FLOW : Approx. Space:")
AdjointFlow_ApproxSpace:print_statistic()

AdjointFlow_ElemDisc = NavierStokesAdjoint("q1,q2,q3, h ", "outer")
AdjointFlow_ElemDisc:set_kinematic_viscosity(visc)
AdjointFlow_ElemDisc:set_stabilization(stab)
AdjointFlow_ElemDisc:set_stabilization_type(stabType)
--Set Imports
--VELOCITY GRADIENT
AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);

--VELOCITY VECTOR
AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)
AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)

--************BOUNDARY CONDITIONS*********--
AdjointFlow_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","inlet")
AdjointFlow_Dirich:add(0,"q2","inlet")
AdjointFlow_Dirich:add(0,"q3","inlet")
--************WALL BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","wall")
AdjointFlow_Dirich:add(0,"q2","wall")
AdjointFlow_Dirich:add(0,"q3","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","obstacle_surface")
AdjointFlow_Dirich:add(0,"q2","obstacle_surface")
AdjointFlow_Dirich:add(0,"q3","obstacle_surface")
--TODO:see if do nothing is ok, might be adjustable think about it...
--************OUTLET BOUNDARY**************--
--AdjointFlow_Dirich:add(0,"q1","outlet")
--AdjointFlow_Dirich:add(0,"q2","outlet")

-- Domain Discretization 
AdjointFlow_DomainDisc = DomainDiscretization(AdjointFlow_ApproxSpace)
AdjointFlow_DomainDisc:add(AdjointFlow_ElemDisc)
AdjointFlow_DomainDisc:add(AdjointFlow_Dirich)

print("ADJOINT FLOW: create GridFunctions and GlobalGridFunctionGradientDatas")
q = GridFunction(AdjointFlow_ApproxSpace);q:set(0.0)
r_q = GridFunction(AdjointFlow_ApproxSpace);r_q:set(0.0)
AdjointFlow_DomainDisc:adjust_solution(q)
q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
q1_value_global=GlobalGridFunctionNumberData(q,"q1")
q2_value_global=GlobalGridFunctionNumberData(q,"q2")
q3_value_global=GlobalGridFunctionNumberData(q,"q3")
---Pressure
h_value_global=GlobalGridFunctionNumberData(q,"h")
---Pressure
A_adjFlow = AssembledLinearOperator(AdjointFlow_DomainDisc)
print("ADJOINT FLOW SYSTEM: all set:")

print("LAGRANGE MATRIX: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
--set up approximation space
Lambda_ApproxSpace = ApproximationSpace(dom)
Lambda_ApproxSpace:add_fct("l1,l2,l3,l4,l5,l6,l7,l8,l9","Piecewise-Constant") 
--Lambda_ApproxSpace:add_fct("l2,l3,l4","Piecewise-Constant") 
Lambda_ApproxSpace:init_levels()
Lambda_ApproxSpace:init_top_surface()

print("LAGRANGE MATRIX: Piecewise constant shape functions Approx. Space:")
Lambda_ApproxSpace:print_statistic()
lambda_piecewise = AdvancedGridFunction(Lambda_ApproxSpace);lambda_piecewise:set(0.0)
q_piecewise = AdvancedGridFunction(Lambda_ApproxSpace);q_piecewise:set(0.0)
q_projected = AdvancedGridFunction(Lambda_ApproxSpace);q_projected:set(0.0)
rhs_piecewise = AdvancedGridFunction(Lambda_ApproxSpace);rhs_piecewise:set(0.0)
temp1_piecewise= AdvancedGridFunction(Lambda_ApproxSpace);temp1_piecewise:set(0.0)--used for difference tau*(Grad_u-q_projected)

lambda00_global=GlobalGridFunctionNumberData(lambda_piecewise,"l1")
lambda01_global=GlobalGridFunctionNumberData(lambda_piecewise,"l2")
lambda02_global=GlobalGridFunctionNumberData(lambda_piecewise,"l3")
lambda10_global=GlobalGridFunctionNumberData(lambda_piecewise,"l4")
lambda11_global=GlobalGridFunctionNumberData(lambda_piecewise,"l5")
lambda12_global=GlobalGridFunctionNumberData(lambda_piecewise,"l6")
lambda20_global=GlobalGridFunctionNumberData(lambda_piecewise,"l7")
lambda21_global=GlobalGridFunctionNumberData(lambda_piecewise,"l8")
lambda22_global=GlobalGridFunctionNumberData(lambda_piecewise,"l9")


--projected value
qproj00_global=GlobalGridFunctionNumberData(q_projected,"l1")
qproj01_global=GlobalGridFunctionNumberData(q_projected,"l2")
qproj02_global=GlobalGridFunctionNumberData(q_projected,"l3")
qproj10_global=GlobalGridFunctionNumberData(q_projected,"l4")
qproj11_global=GlobalGridFunctionNumberData(q_projected,"l5")
qproj12_global=GlobalGridFunctionNumberData(q_projected,"l6")
qproj20_global=GlobalGridFunctionNumberData(q_projected,"l7")
qproj21_global=GlobalGridFunctionNumberData(q_projected,"l8")
qproj22_global=GlobalGridFunctionNumberData(q_projected,"l9")

print("DEFORMATION EQUATION: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
DeformationSpace_ApproxSpace = ApproximationSpace(dom)
DeformationSpace_ApproxSpace:add_fct("u1,u2,u3","Lagrange",1)
DeformationSpace_ApproxSpace:init_levels()
DeformationSpace_ApproxSpace:init_top_surface()
print("DEFORMATION EQUATION: Approx. Space:")
DeformationSpace_ApproxSpace:print_statistic()

--create all necessary grid functions
delta_u = AdvancedGridFunction(DeformationSpace_ApproxSpace);delta_u:set(0.0)
u = AdvancedGridFunction(DeformationSpace_ApproxSpace);u:set(0.0);
u_converged = AdvancedGridFunction(DeformationSpace_ApproxSpace);u_converged:set(0.0);
u_diff = AdvancedGridFunction(DeformationSpace_ApproxSpace);u_diff:set(0.0);
u_old = AdvancedGridFunction(DeformationSpace_ApproxSpace);u_old:set(0.0);
u_negative = AdvancedGridFunction(DeformationSpace_ApproxSpace);u_negative:set(0.0)--to revert transformation
sigma = AdvancedGridFunction(DeformationSpace_ApproxSpace);sigma:set(0.0);
Lu = AdvancedGridFunction(DeformationSpace_ApproxSpace);Lu:set(0.0);
u_zeros = AdvancedGridFunction(DeformationSpace_ApproxSpace);u_zeros:set(0.0)
u1_gradient_global=GlobalGridFunctionGradientData(u,"u1")
u2_gradient_global=GlobalGridFunctionGradientData(u,"u2")
u3_gradient_global=GlobalGridFunctionGradientData(u,"u3")
u1_value_global=GlobalGridFunctionNumberData(u,"u1")
u2_value_global=GlobalGridFunctionNumberData(u,"u2")
u3_value_global=GlobalGridFunctionNumberData(u,"u3")

print("DEFORMATION EQUATION SECOND DERIVATIVE WITHOUT J'' ")
--This is the Linear Operator, aka Hessian, Stiffness Matrix, the lhs
DeformationEquationHessian_ElemDisc = DeformationEquation("u1,u2,u3", "outer");DeformationEquationHessian_ElemDisc:set_quad_order(1)
DeformationEquationHessian_ElemDisc:set_lambda_vol(lambda_vol)
DeformationEquationHessian_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
DeformationEquationHessian_ElemDisc:set_step_length(step_length)
--Set Imports
--DEFORMATION VECTOR
DeformationEquationHessian_ElemDisc:set_deformation_d1(u1_value_global)
DeformationEquationHessian_ElemDisc:set_deformation_d2(u2_value_global)
DeformationEquationHessian_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
DeformationEquationHessian_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
DeformationEquationHessian_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
DeformationEquationHessian_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
--For the linear problem this is the RHS, the Sensitivity gets added later on
DeformationEquationRHS_ElemDisc = DeformationEquationRHS("u1,u2,u3", "outer");DeformationEquationRHS_ElemDisc:set_quad_order(1);
DeformationEquationRHS_ElemDisc:set_lambda_vol(lambda_vol)
DeformationEquationRHS_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
DeformationEquationRHS_ElemDisc:set_step_length(step_length)
DeformationEquationRHS_ElemDisc:set_tau(tau)
--Set Imports
--DEFORMATION VECTOR
DeformationEquationRHS_ElemDisc:set_deformation_d1(u1_value_global)
DeformationEquationRHS_ElemDisc:set_deformation_d2(u2_value_global)
DeformationEquationRHS_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
DeformationEquationRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
DeformationEquationRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
DeformationEquationRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)

--PIECEWISE CONSTANT LAMBDA
DeformationEquationRHS_ElemDisc:set_lambda00(lambda00_global)
DeformationEquationRHS_ElemDisc:set_lambda01(lambda01_global)
DeformationEquationRHS_ElemDisc:set_lambda02(lambda02_global)
DeformationEquationRHS_ElemDisc:set_lambda10(lambda10_global)
DeformationEquationRHS_ElemDisc:set_lambda11(lambda11_global)
DeformationEquationRHS_ElemDisc:set_lambda12(lambda12_global)
DeformationEquationRHS_ElemDisc:set_lambda20(lambda20_global)
DeformationEquationRHS_ElemDisc:set_lambda21(lambda21_global)
DeformationEquationRHS_ElemDisc:set_lambda22(lambda22_global)

--PIECEWISE CONSTANT Q, PROJECTED
DeformationEquationRHS_ElemDisc:set_q00(qproj00_global)
DeformationEquationRHS_ElemDisc:set_q01(qproj01_global)
DeformationEquationRHS_ElemDisc:set_q02(qproj02_global)
DeformationEquationRHS_ElemDisc:set_q10(qproj10_global)
DeformationEquationRHS_ElemDisc:set_q11(qproj11_global)
DeformationEquationRHS_ElemDisc:set_q12(qproj12_global)
DeformationEquationRHS_ElemDisc:set_q20(qproj20_global)
DeformationEquationRHS_ElemDisc:set_q21(qproj21_global)
DeformationEquationRHS_ElemDisc:set_q22(qproj22_global)

--************BOUNDARY CONDITIONS*********--
DeformationEquation_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
DeformationEquation_Dirich:add(0,"u1","inlet")
DeformationEquation_Dirich:add(0,"u2","inlet")
DeformationEquation_Dirich:add(0,"u3","inlet")
--************WALL BOUNDARY**************--
DeformationEquation_Dirich:add(0,"u1","wall")
DeformationEquation_Dirich:add(0,"u2","wall")
DeformationEquation_Dirich:add(0,"u3","wall")
--************OUTLET BOUNDARY**************--
DeformationEquation_Dirich:add(0,"u1","outlet")
DeformationEquation_Dirich:add(0,"u2","outlet")
DeformationEquation_Dirich:add(0,"u3","outlet")

-- Domain Discretization 
DeformationEquation_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
DeformationEquation_DomainDisc:add(DeformationEquationHessian_ElemDisc)--matrix
DeformationEquation_DomainDisc:add(DeformationEquation_Dirich)--boundary conditions
DeformationEquation_DomainDisc:add(DeformationEquationRHS_ElemDisc)--rhs vector without J'
--Set boundary conditions
DeformationEquation_DomainDisc:adjust_solution(sigma)
DeformationEquation_DomainDisc:adjust_solution(u)
A_u_Hessian = AssembledLinearOperator(DeformationEquation_DomainDisc)
print("DEFORMATION EQUATION: all set:")

print("DEFORMATION EQUATION FOR LARGE PROBLEM DISCRETIZATION")
--Create a RHS element discretization
DeformationEquationLargeProblemRHS_ElemDisc = DeformationEquationLargeProblemRHS("u1,u2,u3", "outer");DeformationEquationLargeProblemRHS_ElemDisc:set_quad_order(1);
DeformationEquationLargeProblemRHS_ElemDisc:set_tau(1.0)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_vol(lambda_vol)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
DeformationEquationLargeProblemRHS_ElemDisc:set_step_length(step_length)
--Set Imports
--DEFORMATION VECTOR
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d1(u1_value_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d2(u2_value_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)

--PIECEWISE CONSTANT LAMBDA
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda00(lambda00_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda01(lambda01_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda02(lambda02_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda10(lambda10_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda11(lambda11_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda12(lambda12_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda20(lambda20_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda21(lambda21_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_lambda22(lambda22_global)

--PIECEWISE CONSTANT Q
DeformationEquationLargeProblemRHS_ElemDisc:set_q00(qproj00_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q01(qproj01_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q02(qproj02_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q10(qproj10_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q11(qproj11_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q12(qproj12_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q20(qproj20_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q21(qproj21_global)
DeformationEquationLargeProblemRHS_ElemDisc:set_q22(qproj22_global)

--Large Problem Domain Discretization
DeformationEquationLargeProblem_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
DeformationEquationLargeProblem_DomainDisc:add(DeformationEquationHessian_ElemDisc)
DeformationEquationLargeProblem_DomainDisc:add(DeformationEquationLargeProblemRHS_ElemDisc)
DeformationEquationLargeProblem_DomainDisc:add(DeformationEquation_Dirich)
A_Large = AssembledLinearOperator(DeformationEquationLargeProblem_DomainDisc)

DeformationEquationLargeProblem_DomainDisc:adjust_solution(delta_u)
DeformationEquationLargeProblem_DomainDisc:adjust_solution(u)
print("DEFORMATION EQUATION FOR LARGE PROBLEM: all set")

print("SENSITIVITY J': CREATED")
Jprime_ElemDisc=Sensitivity("u1,u2,u3","outer")
Jprime_ElemDisc:set_kinematic_viscosity(visc)
Jprime_ElemDisc:set_step_length(scaling)
--VELOCITY GRADIENT
Jprime_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
Jprime_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
Jprime_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
--VELOCITY VECTOR
Jprime_ElemDisc:set_velocity_d1(v1_value_global)
Jprime_ElemDisc:set_velocity_d2(v2_value_global)
Jprime_ElemDisc:set_velocity_d3(v3_value_global)
--PRESSURE
Jprime_ElemDisc:set_pressure(p_value_global)
--ADJOINT VELOCITY VECTOR
Jprime_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
Jprime_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
Jprime_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
--ADJOINT VELOCITY GRADIENT
Jprime_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
Jprime_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
Jprime_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
--ADJOINT PRESSURE
Jprime_ElemDisc:set_adjoint_pressure(h_value_global)

Jprime_DomainDisc=DomainDiscretization(DeformationSpace_ApproxSpace)
Jprime_DomainDisc:add(Jprime_ElemDisc)
SensitivityGF=AdvancedGridFunction(DeformationSpace_ApproxSpace)

print("WE MUST CREATE A DISCRETIZATION FOR EACH CONSTRAINT FOR ITS M-SOLVE")
--B-MATRIX VECTORS
Bvol=AdvancedGridFunction(DeformationSpace_ApproxSpace);Bvol:set(0.0);
Bx = AdvancedGridFunction(DeformationSpace_ApproxSpace);Bx:set(0.0);
By = AdvancedGridFunction(DeformationSpace_ApproxSpace);By:set(0.0);
Bz = AdvancedGridFunction(DeformationSpace_ApproxSpace);Bz:set(0.0);

B_vector = {Bvol, Bx, By, Bz}
print("FOR M=1, Bvol")
--Volume derviative
BVolume_ElemDisc = VolumeConstraintSecondDerivative("u1,u2,u3", "outer");BVolume_ElemDisc:set_quad_order(1)
BVolume_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
BVolume_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
BVolume_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
BVolume_ElemDisc:set_deformation_d1(u1_value_global)
BVolume_ElemDisc:set_deformation_d2(u2_value_global)
BVolume_ElemDisc:set_deformation_d3(u3_value_global)

--Bvol domain discretization
Bvol_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
Bvol_DomainDisc:add(BVolume_ElemDisc)
Bvol_DomainDisc:add(DeformationEquationHessian_ElemDisc)
Bvol_DomainDisc:add(DeformationEquation_Dirich)
Avol = AssembledLinearOperator(Bvol_DomainDisc)
t_vol = AdvancedGridFunction(DeformationSpace_ApproxSpace);t_vol:set(0.0);
Bvol_DomainDisc:adjust_solution(t_vol)
print("FOR M=2, Bx")
--Barycenter on X direction
XBarycenter_ElemDisc = SecondDerivativeBarycenter("u1,u2,u3", "outer");XBarycenter_ElemDisc:set_quad_order(1)
XBarycenter_ElemDisc:set_index(1)
XBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
XBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
XBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
XBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
XBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
XBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

--Bx domain discretization
Bx_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
Bx_DomainDisc:add(XBarycenter_ElemDisc)
Bx_DomainDisc:add(DeformationEquationHessian_ElemDisc)
Bx_DomainDisc:add(DeformationEquation_Dirich)
Ax = AssembledLinearOperator(Bx_DomainDisc)
t_x = AdvancedGridFunction(DeformationSpace_ApproxSpace);t_x:set(0.0);
Bx_DomainDisc:adjust_solution(t_x)
print("FOR M=3, Bx")
--Barycenter on Y direction
YBarycenter_ElemDisc = SecondDerivativeBarycenter("u1,u2,u3", "outer");YBarycenter_ElemDisc:set_quad_order(1)
YBarycenter_ElemDisc:set_index(2)
YBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
YBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
YBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
YBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
YBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
YBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

--By domain discretization
By_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
By_DomainDisc:add(YBarycenter_ElemDisc)
By_DomainDisc:add(DeformationEquationHessian_ElemDisc)
By_DomainDisc:add(DeformationEquation_Dirich)
Ay = AssembledLinearOperator(By_DomainDisc)
t_y = AdvancedGridFunction(DeformationSpace_ApproxSpace);t_y:set(0.0);
By_DomainDisc:adjust_solution(t_y)

print("FOR M=4, Bz")
--Barycenter on Z direction
ZBarycenter_ElemDisc = XBarycenterConstraintSecondDerivative("u1,u2,u3", "outer");ZBarycenter_ElemDisc:set_quad_order(1) 
ZBarycenter_ElemDisc:set_index(3)
ZBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
ZBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
ZBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
ZBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
ZBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
ZBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
--By domain discretization
Bz_DomainDisc = DomainDiscretization(DeformationSpace_ApproxSpace)
Bz_DomainDisc:add(ZBarycenter_ElemDisc)
Bz_DomainDisc:add(DeformationEquationHessian_ElemDisc)
Bz_DomainDisc:add(DeformationEquation_Dirich)
Az = AssembledLinearOperator(Bz_DomainDisc)
t_z = AdvancedGridFunction(DeformationSpace_ApproxSpace);t_z:set(0.0);
Bz_DomainDisc:adjust_solution(t_z)

print("CREATION OF DATA STRUCTURES FOR THE SCHUR COMPLEMENT, BASED ON MATRIX PACKAGE")
BTranspose_sigma = matrix(m,1)
L_lambda = matrix(m,1)--format from package lua-matrix
Lambda = matrix(m,1)
Lambda[1][1]=lambda_vol
Lambda[2][1]=lambda_x
Lambda[3][1]=lambda_y
Lambda[4][1]=lambda_z

S = matrix(m,m)
inverse_S = matrix(m,m)
rhs = matrix(m,1)
DeltaLambda=matrix(m,1)
NegativeB = AdvancedGridFunction(DeformationSpace_ApproxSpace);NegativeB:set(0.0);
MinusLu_BdeltaLambda = AdvancedGridFunction(DeformationSpace_ApproxSpace);MinusLu_BdeltaLambda:set(0.0);


print("CREATION OF MASS MODEL WITH A DIAG(MASS MATRTIX")
print("The appropriate vectors are created before this point")
MassModel_ElemDisc = MassModel("l1,l2,l3,l4,l5,l6,l7,l8,l9","outer")
--Set Imports
--DEFORMATION VECTOR
MassModel_ElemDisc:set_deformation_d1(u1_value_global)
MassModel_ElemDisc:set_deformation_d2(u2_value_global)
MassModel_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
MassModel_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
MassModel_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
MassModel_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
MassModel_ElemDisc:set_lambda00(lambda00_global)
MassModel_ElemDisc:set_lambda01(lambda01_global)
MassModel_ElemDisc:set_lambda02(lambda02_global)
MassModel_ElemDisc:set_lambda10(lambda10_global)
MassModel_ElemDisc:set_lambda11(lambda11_global)
MassModel_ElemDisc:set_lambda12(lambda12_global)
MassModel_ElemDisc:set_lambda20(lambda20_global)
MassModel_ElemDisc:set_lambda21(lambda21_global)
MassModel_ElemDisc:set_lambda22(lambda22_global)
-- Domain Discretization 
MassModel_DomainDisc = DomainDiscretization(Lambda_ApproxSpace)
MassModel_DomainDisc:add(MassModel_ElemDisc)
DiagQ = AssembledLinearOperator(MassModel_DomainDisc)

--This discretization performs tau*(Grad_u-q_projected)
LambdaUpdate_ElemDisc=LambdaUpdate("l1,l2,l3,l4,l5,l6,l7,l8,l9","outer")
--Set Imports
--DEFORMATION GRADIENT
LambdaUpdate_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
LambdaUpdate_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
LambdaUpdate_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
LambdaUpdate_ElemDisc:set_qproj00(qproj00_global)
LambdaUpdate_ElemDisc:set_qproj01(qproj01_global)
LambdaUpdate_ElemDisc:set_qproj02(qproj02_global)
LambdaUpdate_ElemDisc:set_qproj10(qproj10_global)
LambdaUpdate_ElemDisc:set_qproj11(qproj11_global)
LambdaUpdate_ElemDisc:set_qproj12(qproj12_global)
LambdaUpdate_ElemDisc:set_qproj20(qproj20_global)
LambdaUpdate_ElemDisc:set_qproj21(qproj21_global)
LambdaUpdate_ElemDisc:set_qproj22(qproj22_global)
-- Domain Discretization 
LambdaUpdate_DomainDisc = DomainDiscretization(Lambda_ApproxSpace)
LambdaUpdate_DomainDisc:add(LambdaUpdate_ElemDisc)


print("CREATION OF SOLVERS")
NavierStokes_Solver=util.oo.ns_solver(NavierStokes_DomainDisc,NavierStokes_ApproxSpace)
AdjointFlow_Solver = util.oo.adjoint_ns_solver(AdjointFlow_DomainDisc,AdjointFlow_ApproxSpace)

ADMMDiagonal_Solver=CG();
ADMMDiagonal_Solver:set_preconditioner(Jacobi(0.66));--without this, we have plenty more iterations
ADMMDiagonal_Solver:set_convergence_check(ConvCheck(2000,1.0e-9,0.0,true))
SmallProblemRHS_Solver = util.oo.linear_solver(DeformationEquation_DomainDisc,DeformationSpace_ApproxSpace,false)
Bvol_Solver = util.oo.linear_solver(Bvol_DomainDisc,DeformationSpace_ApproxSpace,false)
Bx_Solver = util.oo.linear_solver(Bx_DomainDisc,DeformationSpace_ApproxSpace,false)
By_Solver = util.oo.linear_solver(By_DomainDisc,DeformationSpace_ApproxSpace,false)
Bz_Solver = util.oo.linear_solver(Bz_DomainDisc,DeformationSpace_ApproxSpace,false)
LargeProblem_Solver = util.oo.linear_solver(DeformationEquationLargeProblem_DomainDisc,DeformationSpace_ApproxSpace,true)
print("NAVIER STOKES SOLVER IS:")
print(NavierStokes_Solver:config_string())
print("ADJOINT FLOWS SOLVER IS:")
print(AdjointFlow_Solver:config_string())

--This writer is used for all vtk file printings     
vtkWriter = VTKOutput()

--Storage of data
vStep = {}
vDrag = {}
vShapeDerivative = {}
vSupNorm = {}
--vStepSize = {}
--vStepLength = {}
--Solver data storage per step
vTotalLinearIterations = {}
vLargeSolverIterations = {}
vBvolSolverIterations = {}
vBxSolverIterations = {}
vBySolverIterations = {}
vBzSolverIterations = {}
vRHSSolver = {}
vNonLinearIterationsPerP = {}

sum_total_iterations = 0
sum_largesolver_iterations = 0
sum_bvolsolver_iterations = 0
sum_bxsolver_iterations = 0
sum_bysolver_iterations = 0
sum_bzsolver_iterations = 0
sum_rhssolver_iterations = 0
sum_newtonsteps = 0
--Lambda multipliers
vLambdaVol ={}
vLambdaX ={}
vLambdaY ={}
vLambdaZ ={}
--variables used across the admm implementation
maximum_frobenius_norm=0.0;
admm_steps=0;
u_diff_norm=0.0;
lambda_inc_norm=0.0

--The first Navier Stokes problem solve is done outside the loop
print("SOLVE PHASE: NON-LINEAR SOLUTION OF THE NAVIER STOKES PROBLEM")
NavierStokes_Solver:init(navier_Op)----TODO:navierLine
if NavierStokes_Solver:prepare(v) == false then print ("NavierStokes: Newton solver prepare failed at step "..step.."."); PrintStats(); exit(); end 			
if NavierStokes_Solver:apply(v) == false then print ("NavierStokes: Newton solver failed at step "..step.."."); PrintStats(); exit(); end 

VecScaleAssign(v_temp,1.0,v)--store as initial guess for next solves, v is used for imports
if bOutputFlows and bReadInitialGuess == -1 then
	print("Flows are saved to '" .. flowOutputFile .. "'...")
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(flowNames, nodalFlows)
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print(flowOutputFile, v, 0,0, false)
end
if bOutputPressure and bReadInitialGuess == -1 then
	--PRESSURE
	print("Pressure field is saved")
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("p", "pressure")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("pressureOutputFile", v,0,0, false)
end

step = 0;p_solver_failure = false
drag_old = 0.5 * visc * Drag(u_zeros, v, "v1,v2,v3","u1,u2,u3","outer",3)
vDrag[step] = 0.5 * visc * Drag(u_zeros, v, "v1,v2,v3","u1,u2,u3","outer",3)
ReferenceVolume = VolumeDefect(u,0, "outer", "u1,u2,u3",4,false,1,false)--used as a constant
vBarycenterDefect = {0,0,0}
print("The reference volume is: "..ReferenceVolume)


--Begin of optimization loop
print("+++++++++++++++++++ BEGIN OPTIMIZATION LOOP  +++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
while (step <= numSteps) do
	vStep[step] = step
	print("		+++++++++++++++++++ NEW OPTIMIZATION STEP "..step.." +++++++++++++++")
	print("		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print("		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

	if bDebugOutput	then SaveGridLevelToFile(dom:grid(), dom:subset_handler(), numRefs, "Mesh_lev"..numRefs.."_step"..step..".ugx") end
	
	print("		SOLVE PHASE: LINEAR SOLUTION OF THE ADJOINT FLOWS PROBLEM")
	AdjointFlow_DomainDisc:assemble_linear(A_adjFlow, r_q)--TODO:adjointLine
	AdjointFlow_Solver:init(A_adjFlow, q)
	if AdjointFlow_Solver:apply(q, r_q) == false then print("Adjoint flows solver failed, at step "..ns_i.." p "..p_current);exit(); end
	if bOutputAdjoints then
		print("VTK_WRITER::Adjoints is saved to '" .. adjointFlowOutputFile .. "'...")
		vtkWriter:clear_selection()
		vtkWriter:select_nodal(adjointFlowNames, adjointNodalFlows)
		--vtkWriter:select_all(false)	-- write all variables
		vtkWriter:print(adjointFlowOutputFile, q,step,step, false)
		--PRESSURE
		vtkWriter:clear_selection()
		vtkWriter:select_nodal("h", "adjoint_pressure")
		--vtkWriter:select_all(false)	-- write all variables
		vtkWriter:print("adjointPressureOutputFile", q,step,step, false)
	end

	--Assemble J'as a standalone vector, it has a special setting of gauss points
	--and its the only vector that gets SetZeroAwayFromSubset'ed
	scaling=util.GetParamNumber("-scaling",1.0);Jprime_ElemDisc:set_step_length(scaling);Jprime_DomainDisc:assemble_defect(SensitivityGF,u)
	SetZeroAwayFromSubset(SensitivityGF,"u1,u2,u3","obstacle_surface")
	--SensitivityGF:change_storage_type_to_consistent()
	if bDebugSensitivity then
		SensitivityGF:change_storage_type_to_consistent()
		print("SensitivityGF is printed")
		vtkWriter:clear_selection()
		vtkWriter:select_nodal("u1,u2,u3", "J'")
		--vtkWriter:select_all(false)	-- write all variables
		vtkWriter:print("sensitivity", SensitivityGF, step,step, false)
		SensitivityGF:change_storage_type_to_additive()
	end
	--VecScaleAssign(SensitivityGF,-1.0,SensitivityGF)

	--variables store totals for each step, set to zero upon starting a new step
	sum_total_iterations = 0
	sum_largesolver_iterations = 0
	sum_bvolsolver_iterations = 0
	sum_bxsolver_iterations = 0
	sum_bysolver_iterations = 0
    sum_bzsolver_iterations = 0
	sum_rhssolver_iterations = 0
	sum_newtonsteps = 0
	

	while(true) do--step size control loop

		--Here we set to zero certain data-structures for the beginning of a new optimization step
		--or in case a step is repeated
		--L_lambda is zero at start of new step
		L_lambda[1][1] = 0.0
		L_lambda[2][1] = 0.0
		L_lambda[3][1] = 0.0
        L_lambda[4][1] = 0.0
		if p_solver_failure==false then
			VecScaleAssign(u_converged,1.0,u);
		end
		u:set(0.0);lambda_piecewise:set(0.0)

		sigma:set(0.0)
		Lambda[1][1]=0.0;Lambda[2][1]=0.0;Lambda[3][1]=0.0;Lambda[4][1]=0.0;
	
		--TODO:when would this be necessary, how does failure look like? In the end, what is failure?
		if p_solver_failure == true then--step had failed, set bool back to to false
			p_solver_failure=false--restart boolean
			print("+++++++++++++++++++ AT OPTIMIZATION STEP "..step.." +++++++++++++++++")
			print("++++++++++++++ THE STEP IS REPEATED BECAUSE OF SOLVER FAILURE+++++++")
			print("+++++++++++++++THE ADJOINTS ARE NOT CALCULATED +++++++++++++++++++++")
		end
		--After this loop we must obtain the deformation field
		--Vectors for information
		vADMM_Step = {}
		vADMM_Scaling = {}
		vADMM_Sigma = {}
		vADMM_Udiff = {}
		vADMM_LambdaInc = {}
		vADMM_SigmaMinusMaxNorm = {}
		vADMM_MaxFrobNorm = {}
		admm_steps=0.0
		while (true and admm_steps<admmSteps) do--ADMM loop
			print("******************START ADMM STEP "..admm_steps.." ******************")
			print("+++++++++++++++OPTIMIZATION STEP "..step.."+++++++++++++++")
			--Some marvelos stuff is supposed to happen here
			--Marvel, wonder and believe in ADMM so that you may have life
			
			--Prepare a new ADMM step
			--u:set(0.0);
			--Lambda[1][1]=0.0;Lambda[2][1]=0.0;Lambda[3][1]=0.0;
			L_lambda[1][1] = 0.0
		    L_lambda[2][1] = 0.0
		    L_lambda[3][1] = 0.0
            L_lambda[4][1] = 0.0

			DeformationEquationHessian_ElemDisc:set_lambda_vol(0.0)
            DeformationEquationHessian_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)
            DeformationEquationRHS_ElemDisc:set_lambda_vol(0.0)
            DeformationEquationRHS_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)
            DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_vol(0.0)
            DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)

            --Mass model solution
            q_piecewise:set(0.0);rhs_piecewise:set(0.0);temp1_piecewise:set(0.0)
			--solve mass system
            MassModel_DomainDisc:assemble_jacobian(DiagQ,u_negative)
            MassModel_DomainDisc:assemble_defect(rhs_piecewise,u_negative)
            --VecScaleAssign(rhs_piecewise,-1.0,rhs_piecewise)
			--MassModel_DomainDisc:assemble_linear(DiagQ, rhs_piecewise)
			ADMMDiagonal_Solver:init(DiagQ, q_piecewise)
			if ADMMDiagonal_Solver:apply(q_piecewise, rhs_piecewise) == false then print("Mass model solver failed at admm step "..admm_steps);p_solver_failure = true;break; end
			VecScaleAssign(q_piecewise,-1.0,q_piecewise)
            --WriteGridFunctionToVTK(q_piecewise, "q_step_"..admm_steps)
            --WriteGridFunctionToVTK(rhs_piecewise, "mass_model_rhs_"..admm_steps)
            --q_projected:set(0.0)
			--threshold q_piecewise
			Testing(q_projected, q_piecewise,"l1,l2,l3,l4,l5,l6,l7,l8,l9",sigma_threshold);
            --exit()
			q_projected:change_storage_type_to_consistent()
			--WriteGridFunctionToVTK(q_projected, "q_projected_step_"..admm_steps)

            --Threshold
            maximum_frobenius_norm=MaximumFrobeniusNorm(u_old,"u1,u2,u3","outer",4)
            print("Maximum Frobenius norm is: "..maximum_frobenius_norm)
		

			--START NEWTON SOLVER METHOD, a simplification of algorithm 2 in [1] for constant p=2
			ns_i=1
			--Store convergence of Newton's method and residual for the current step
			vNS_Step = {}
			vNS_NormDeltaUp = {}
			vNS_NormDeltaLambda= {}
			vNS_NormSum = {}
			vNS_Lu2Norm = {}
			
			vNS_RHS_Solver_Iterations = {}
			vNS_LargeSolver_Iterations = {}
			vNS_BxSolver_Iterations = {}
			vNS_BySolver_Iterations = {}
            vNS_BzSolver_Iterations = {}
			vNS_BvolSolver_Iterations = {}

			--print("******************START: NEWTON'S METHOD******************")
			--print("+++++++++++++++OPTIMIZATION STEP "..step.."+++++++++++++++")
			ns_i = 1--set newton step counter

			while(ns_i <= nsMaxIts) do--Newton's Method for nonlinear system of equations
				--Data structures used for Newton's method convergence-related files
				vNS_Step[ns_i]=ns_i
				--Data structures are set to 0.0 at each new nonlinear iteration
				MinusLu_BdeltaLambda:set(0.0);Lu:set(0.0)
                delta_u:set(0.0);sigma:set(0.0);
				rhs[1][1]=0.0;rhs[2][1]=0.0;rhs[3][1]=0.0;rhs[4][1]=0.0;
				DeltaLambda[1][1]=0.0;DeltaLambda[2][1]=0.0;DeltaLambda[3][1]=0.0;DeltaLambda[4][1]=0.0;
				BTranspose_sigma[1][1]=0.0;BTranspose_sigma[2][1]=0.0;BTranspose_sigma[3][1]=0.0;BTranspose_sigma[4][1]=0.0;
				--SensitivityGF:set(0.0);Jprime_DomainDisc:assemble_defect(SensitivityGF,u)
				--forces update of u in the assemby
				--at set(0.0) the gfs are PST_CONSISTENT
				Bvol:set(0.0);Bx:set(0.0);By:set(0.0);Bz:set(0.0);
				--after assemble_defect the gfs are PST_ADDITIVE
				Bvol_DomainDisc:assemble_defect(Bvol,u);Bx_DomainDisc:assemble_defect(Bx,u);By_DomainDisc:assemble_defect(By,u);Bz_DomainDisc:assemble_defect(Bz,u)
				--these correct the minus sign introduced by the assemble_defect routine
				VecScaleAssign(Bvol,-1.0,Bvol)
				VecScaleAssign(Bx,-1.0,Bx)
				VecScaleAssign(By,-1.0,By)
                VecScaleAssign(Bz,-1.0,Bz)

				B_vector = {Bvol, Bx, By, Bz}
				print("			****************** NEWTON'S METHOD ITERATION 	#:"..ns_i.." *****************************")
				print("			****************** OPTIMIZATION STEP:"..step.." ******************************************")
				print("			****************** ADMM STEP:"..admm_steps.." ******************************************")
				print("			****************** SCALING VALUE IS: "..scaling.." *******************")
				
				--We start by solving the second equation from (28) in [1]
				--corresponding to algorithm 3 in [1]
				--Calculate rhs = -L_lambda+BTranspose.sigma, solve A.sigma=Lu for sigma
				print("			SCHUR PROBLEM:: 1.- SOLVE LINEAR PROBLEM OF RHS: A.sigma=(-Lu)")	
				DeformationEquation_DomainDisc:adjust_solution(sigma)--PST_CONSISTENT
                DeformationEquation_DomainDisc:assemble_jacobian(A_u_Hessian, u)
                DeformationEquation_DomainDisc:assemble_defect(Lu,u)--Lu is additive
				--VecScaleAssign(Lu,1.0,Lu)
				--SetZeroAwayFromSubset(Lu,"u1,u2","obstacle_surface")
				VecScaleAdd2(Lu,1.0,Lu,1.0,SensitivityGF)

				if Lu:has_storage_type_additive() == false then print("CATASTROPHIC FAILURE::RHS NOT ADDITIVE");exit();end
				SmallProblemRHS_Solver:init(A_u_Hessian, sigma)
				if SmallProblemRHS_Solver:apply(sigma, Lu) == false then print("A.sigma=Lu, solver failed, at step "..ns_i);p_solver_failure = true; break; end
                VecScaleAssign(sigma,-1.0,sigma)
				Lu:change_storage_type_to_consistent()--this ruins Lu
				if bDebugOutput then
					print("Lu is saved to ConsistentLu_step_"..step)
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2","up")
					vtkWriter:print("ConsistentLu_step_"..step, Lu,ns_i,ns_i, false)
				end

				--Perform B.delta_lambda
				--sigma is consistent while B_vector[i] is additive, this is ok
				--in vector prod <left, right> left changes of storage if necessary
				for i = 1, m do
					BTranspose_sigma[i][1]=VecProd(B_vector[i], sigma);--This is B.sigma
				end
				--assembly of rhs for small problem
				--this step control should do nothing
				rhs[1][1]= -1.0*L_lambda[1][1] - BTranspose_sigma[1][1]
				rhs[2][1]= -1.0*L_lambda[2][1] - BTranspose_sigma[2][1]
				rhs[3][1]= -1.0*L_lambda[3][1] - BTranspose_sigma[3][1]
                rhs[4][1]= -1.0*L_lambda[4][1] - BTranspose_sigma[4][1]

				rhs:print()
				--Volume Constraint...set(0.0) makes PST_CONSISTENT, but ass_def makes ADDITIVE
				Bvol:set(0.0)
                Bvol_DomainDisc:assemble_defect(Bvol,u)
				--SetZeroAwayFromSubset(Bvol,"u1,u2","obstacle_surface")
				Bvol_DomainDisc:assemble_jacobian(Avol,u)
				Bvol_Solver:init(Avol, t_vol)
				if Bvol_Solver:apply(t_vol, Bvol) == false then print("Solver for B[1] failed");p_solver_failure = true;break; end
				--VecScaleAssign(t_vol,-1.0,t_vol)
				VecScaleAssign(Bvol,-1.0,Bvol)
	
				S[1][1]= VecProd(Bvol,t_vol)
				S[2][1]= VecProd(Bx,t_vol)
				S[3][1]= VecProd(By,t_vol)
                S[4][1]= VecProd(Bz,t_vol)

				
				--X-Barycenter Constraint
				Bx:set(0.0)
                Bx_DomainDisc:assemble_defect(Bx,u)
				--SetZeroAwayFromSubset(Bx,"u1,u2","obstacle_surface")
				Bx_DomainDisc:assemble_jacobian(Ax,u)
				Bx_Solver:init(Ax, t_x)
				if Bx_Solver:apply(t_x, Bx) == false then print("Solver for B[2] failed");p_solver_failure = true; break; end	
				-- VecScaleAssign(t_x,-1.0,t_x)
				VecScaleAssign(Bx,-1.0,Bx)

				S[1][2]= VecProd(Bvol,t_x)
				S[2][2]= VecProd(Bx,t_x)
				S[3][2]= VecProd(By,t_x)
                S[4][2]= VecProd(Bz,t_x)
			
				--Y-Barycenter Constraint
				By:set(0.0)
				By_DomainDisc:assemble_defect(By,u)
				--SetZeroAwayFromSubset(By,"u1,u2","obstacle_surface")
				By_DomainDisc:assemble_jacobian(Ay,u)
				By_Solver:init(Ay, t_y)
				if By_Solver:apply(t_y, By) == false then print("Solver for B[3] failed");p_solver_failure = true; break; end	
				VecScaleAssign(By,-1.0,By)

				S[1][3]= VecProd(Bvol,t_y)
				S[2][3]= VecProd(Bx,t_y)
				S[3][3]= VecProd(By,t_y)
                S[4][3]= VecProd(Bz,t_y)

                --Z-Barycenter Constraint
				Bz:set(0.0)
                Bz_DomainDisc:assemble_defect(Bz,u)
				Bz_DomainDisc:assemble_jacobian(Az,u)
				Bz_Solver:init(Az, t_z)
				if Bz_Solver:apply(t_z, Bz) == false then print("Solver for B[3] failed");p_solver_failure = true; break; end	
                VecScaleAssign(Bz,-1.0,Bz)
				S[1][4]= VecProd(Bvol,t_z)
				S[2][4]= VecProd(Bx,t_z)
				S[3][4]= VecProd(By,t_z)
				S[4][4]= VecProd(Bz,t_z)

				print("		THE SCHUR COMPLEMENT MATRIX S IS:")	
				S:print()
				inverse_S = S:invert()
				print("		THE INVERSE OF SCHUR COMPLEMENT MATRIX S IS:")	
				inverse_S:print()
				print("		<S,inv(S)> RETURNS IDENTITY MATRIX:")	
				Iden = matrix(m,m);Iden = S:mul(inverse_S);Iden:print();
		
				--solve reduced problem by direct inverse
				DeltaLambda=matrix(m,1)
				DeltaLambda=inverse_S:mul(rhs)
				print("DeltaLambda IS:")
				DeltaLambda:print()

				lhs = matrix(m,1)
				lhs= S:mul(DeltaLambda)
				diff = matrix(m,1)
				diff=lhs:sub(rhs)
				
				--solve big problem	
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_vol(DeltaLambda[1][1])
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_bx(DeltaLambda[2][1])
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_by(DeltaLambda[3][1])
                DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_bz(DeltaLambda[4][1])

				MinusLu_BdeltaLambda:set(0.0)--does setting it to 0.0 correct the storage type?
				DeformationEquationLargeProblem_DomainDisc:adjust_solution(delta_u)
				print("		LARGE PROBLEM SOLVER AT OPTIMIZATION STEP "..step.." AND NS STEP#   "..ns_i.." ")

				DeformationEquationLargeProblem_DomainDisc:assemble_jacobian(A_Large, u_zeros)--pLaplaceLine
				DeformationEquationLargeProblem_DomainDisc:assemble_defect(MinusLu_BdeltaLambda,u_zeros)
				VecScaleAdd2(MinusLu_BdeltaLambda,1.0,MinusLu_BdeltaLambda,1.0,SensitivityGF)

				LargeProblem_Solver:init(A_Large, delta_u)
				if LargeProblem_Solver:apply_return_defect(delta_u, MinusLu_BdeltaLambda) == false then print("Large problem solver failed at step "..ns_i);p_solver_failure = true; PrintStats();break; end		
				MinusLu_BdeltaLambda:change_storage_type_to_consistent()--this ruins the vector
				if bDebugOutput then
					print("-Lu-B.delta_lambda is saved to RHSBigProb_"..step.."")
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2,u3","up")
					vtkWriter:print("RHSBigProb_"..step.."", MinusLu_BdeltaLambda,ns_i,ns_i, false)
				end
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_vol(0.0)
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_bx(0.0)
				DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_by(0.0)
                DeformationEquationLargeProblemRHS_ElemDisc:set_multiplier_bz(0.0)
				DeformationEquationLargeProblem_DomainDisc:adjust_solution(delta_u)
				--UPDATE DEFORMATION ITERATE u_new=u_old+delta_u
				VecScaleAdd2(u,1.0,u,-1.0,delta_u)
				if bDebugOutput then
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2","u")
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print("delta_u", delta_u,ns_i,ns_i, false)
				end
				--UPDATE LAGR.MULTS. ITERATE lambda_new=lambda_old+lambda_u
				Lambda[1][1]= Lambda[1][1] + DeltaLambda[1][1]
				Lambda[2][1]= Lambda[2][1] + DeltaLambda[2][1] 
				Lambda[3][1]= Lambda[3][1] + DeltaLambda[3][1]	
                Lambda[4][1]= Lambda[4][1] + DeltaLambda[4][1]	
				--re adjust u_p
				DeformationEquation_DomainDisc:adjust_solution(u)
				--counter increase and convergence check
				ns_i = ns_i+1

				if ns_i > nsMaxIts then
					print("		++++++++++++++++ NEWTON METHOD DID NOT CONVERGE:  ++++++++++++++++++++++++++")
					print("		++++++++++++++++ MAX NUMBER OF ITERATIONS REACHED ++++++++++++++++++++++++++")
					print("		++++++++++++++++ AT ITERATION   #: "..ns_i.." ++++++++++++++++++++++++++++")
					print("		++++++++++++++++ AND OPT. STEP # : "..step.." ++++++++++++++++++++++++++++")
					p_solver_failure=true;break;--breaks from newton loop into p loop
				end
				
				--norm of the residual Lu
				--lu_norm1=VecNorm(Lu)
				--lu_norm2=0.0
				lu_norm1=L2Norm(Lu,"u1",4,"outer")
				lu_norm2=L2Norm(Lu,"u2",4,"outer")
                lu_norm3=L2Norm(Lu,"u3",4,"outer")
				lu_norm_sum = math.sqrt(lu_norm1*lu_norm1+lu_norm2*lu_norm2+lu_norm3*lu_norm3)
			
				--norm of the 0th delta_up iterate, starts at 0, could also be overriden
				delta_u_1=L2Norm(delta_u,"u1",4,"outer")
				delta_u_2=L2Norm(delta_u,"u2",4,"outer")
                delta_u_3=L2Norm(delta_u,"u3",4,"outer")
				delta_u_norm_sum=math.sqrt(delta_u_1*delta_u_1+delta_u_2*delta_u_2+delta_u_3*delta_u_3)
				--euclidean norm of delta_lambda, starts at 0
				delta_lambda_norm=0.0
				for j=1,m do
					delta_lambda_norm= delta_lambda_norm+DeltaLambda[j][1]*DeltaLambda[j][1]
				end
				delta_lambda_norm = math.sqrt(delta_lambda_norm)
				delta_norms_sum= delta_lambda_norm--+delta_u_p_norm_sum

				vNS_NormDeltaUp[ns_i-1] = delta_u_norm_sum
				vNS_NormDeltaLambda[ns_i-1] = delta_lambda_norm
				vNS_NormSum[ns_i-1] = 0.0--delta_norms_sum
				vNS_Lu2Norm[ns_i-1]=lu_norm_sum;
				--for iteration counts for every p in every step
				vNS_RHS_Solver_Iterations[ns_i-1]= SmallProblemRHS_Solver:step();
				vNS_LargeSolver_Iterations[ns_i-1]= LargeProblem_Solver:step()
				vNS_BxSolver_Iterations[ns_i-1] = Bx_Solver:step()
				vNS_BySolver_Iterations[ns_i-1] = By_Solver:step()
                vNS_BzSolver_Iterations[ns_i-1] = Bz_Solver:step()
				vNS_BvolSolver_Iterations[ns_i-1] = Bvol_Solver:step()
				--compute the geometrical constraints, assign to vector, negative is used above
				L_lambda[1][1]=VolumeDefect(u, ReferenceVolume, "outer", "u1,u2,u3",4,false,1,false)
				vBarycenterDefect = BarycenterDefect(u,"u1,u2,u3","outer",4)
				L_lambda[2][1]= vBarycenterDefect[1]
				L_lambda[3][1]= vBarycenterDefect[2]
                L_lambda[4][1]= vBarycenterDefect[3]
				print("		GEOMETRICAL CONSTRAINTS DEFECTS (L_LAMBDA) AT NS STEP#   "..ns_i.."  ARE:")	
				L_lambda:print()
				
				--prepare for next iteration
				--Update Lagrange multipliers of the RHS of A.sigma=Lu, all A matrices, and RHS of Large problem
				--Hessian Matrix
				DeformationEquationHessian_ElemDisc:set_lambda_vol(Lambda[1][1])
				DeformationEquationHessian_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				--RHS of Small Problem
				DeformationEquationRHS_ElemDisc:set_lambda_vol(Lambda[1][1 ])
				DeformationEquationRHS_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				--RHS of Large Problem
				DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_vol(Lambda[1][1])
				DeformationEquationLargeProblemRHS_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				
				--iteration counts update
				sum_rhssolver_iterations = sum_rhssolver_iterations + SmallProblemRHS_Solver:step()
				sum_largesolver_iterations = sum_largesolver_iterations + LargeProblem_Solver:step()
				sum_bvolsolver_iterations = sum_bvolsolver_iterations + Bvol_Solver:step()
				sum_bxsolver_iterations = sum_bxsolver_iterations + Bx_Solver:step()
				sum_bysolver_iterations = sum_bysolver_iterations + By_Solver:step()
                sum_bzsolver_iterations = sum_bzsolver_iterations + Bz_Solver:step()


				print("#			"..ns_i.." INCREMENT NORM IS:    "..delta_norms_sum)
				print("#			"..ns_i.." DELTA_U INCREMENT NORM IS:    "..delta_u_norm_sum.."   DELTA_LAMBDA INCREMENT NORM IS:   "..delta_lambda_norm)
				if delta_norms_sum <= nsTol then
					break;--break from newton into p loop
				end

			end--end Newton's Method while loop

			--check immediately convergence and solver failure
			if(p_solver_failure==true) then 
				print("ADMM LOOP::solver failure, breaking");
				break;
			end
			--[[
			if true then
				vtkWriter:clear_selection()
				vtkWriter:select_nodal("u1,u2","u")
				--vtkWriter:select_all(false)	-- write all variables
				vtkWriter:print("u_admmStep", u,admm_steps,admm_steps, false)
			end
			]]
			--exit()

			--Calculate lambda
			--lambda += tau(Grad_u-q_proj)
			LambdaUpdate_DomainDisc:assemble_defect(temp1_piecewise,u_negative);VecScaleAssign(temp1_piecewise, -1.0, temp1_piecewise)
			temp1_piecewise:change_storage_type_to_consistent()
			VecScaleAdd2(lambda_piecewise,1.0,lambda_piecewise,1.0,temp1_piecewise);
			--lambda_piecewise:set(0.0)
			lambda_piecewise:change_storage_type_to_consistent()
            --WriteGridFunctionToVTK(lambda_piecewise, "lambda_step_"..admm_steps)
			--FlattenGradient(q,lambda_piecewise,"l1,l2,l3,l4",u,"u1,u2","outer")
			--exit()
			--u has the new solution
			--u_old has the old one
            u_diff:set(0.0)
			VecScaleAdd2(u_diff,1.0,u,-1.0,u_old);--this is used as a convergence measure
			VecScaleAssign(u_old,1.0,u);
            

			--start convergence check preparations
			u_diff_norm1=L2Norm(u_diff,"u1",4,"outer")
			u_diff_norm2=L2Norm(u_diff,"u2",4,"outer")
            u_diff_norm3=L2Norm(u_diff,"u3",4,"outer")
			u_diff_norm = math.sqrt(u_diff_norm1*u_diff_norm1+u_diff_norm2*u_diff_norm2+u_diff_norm3*u_diff_norm3)
			--u_diff_norm=VecNorm(u_diff)
			--TODO:l2norm for increment
			lambda_diff_norm1=L2Norm(temp1_piecewise,"l1",4,"outer")
			lambda_diff_norm2=L2Norm(temp1_piecewise,"l2",4,"outer")
			lambda_diff_norm3=L2Norm(temp1_piecewise,"l3",4,"outer")
			lambda_diff_norm4=L2Norm(temp1_piecewise,"l4",4,"outer")
            lambda_diff_norm5=L2Norm(temp1_piecewise,"l5",4,"outer")
            lambda_diff_norm6=L2Norm(temp1_piecewise,"l6",4,"outer")
            lambda_diff_norm7=L2Norm(temp1_piecewise,"l7",4,"outer")
            lambda_diff_norm8=L2Norm(temp1_piecewise,"l8",4,"outer")
            lambda_diff_norm9=L2Norm(temp1_piecewise,"l9",4,"outer")
          
			lambda_inc_norm=math.sqrt(lambda_diff_norm1*lambda_diff_norm1+lambda_diff_norm2*lambda_diff_norm2+lambda_diff_norm3*lambda_diff_norm3+lambda_diff_norm4*lambda_diff_norm4+lambda_diff_norm5*lambda_diff_norm5+lambda_diff_norm6*lambda_diff_norm6+lambda_diff_norm7*lambda_diff_norm7+lambda_diff_norm8*lambda_diff_norm8+lambda_diff_norm9*lambda_diff_norm9)
			--lambda_inc_norm=VecNorm(temp1_piecewise)

			print("ADMM LOOP::STEP="..admm_steps)
			print("ADMM loop convergence check: MaxFrobNorm= "..maximum_frobenius_norm)
			print("ADMM loop convergence check: sigma-max(Grad(u))= "..sigma_threshold-maximum_frobenius_norm)
			print("ADMM loop convergence check: abs(max(Grad(u))-sigma)= "..math.abs(maximum_frobenius_norm-sigma_threshold))
			print("ADMM loop convergence check: L2Norm of u_diff= "..u_diff_norm)
			print("ADMM loop convergence check: VecNorm of lambda_inc= ".. lambda_inc_norm)
			print("ADMM loop convergence check: scaling= "..scaling)
			print("ADMM loop convergence check: threshold= "..sigma_threshold)

			vADMM_Step[admm_steps]=admm_steps
			--table.insert(vADMM_Step,admm_steps)
			vADMM_Scaling[admm_steps]=scaling
			vADMM_Sigma[admm_steps]=sigma_threshold
			vADMM_Udiff[admm_steps]=u_diff_norm
			vADMM_LambdaInc[admm_steps]=lambda_inc_norm
			vADMM_MaxFrobNorm[admm_steps]=maximum_frobenius_norm
			vADMM_SigmaMinusMaxNorm[admm_steps]=sigma_threshold-maximum_frobenius_norm
			 
			gnuplot.write_data("__ADMMStats_step_"..step.."_.txt", {vADMM_Step,vADMM_Scaling,vADMM_Sigma,
																		vADMM_Udiff,vADMM_LambdaInc,
																		vADMM_MaxFrobNorm,vADMM_SigmaMinusMaxNorm},false)

			--CONVERGENCE CHECK OF ADMM LOOP
			if((lambda_inc_norm < admm_tolerance) and (u_diff_norm< admm_tolerance) and (sigma_threshold - maximum_frobenius_norm > -0.05*sigma_threshold)) then --check if converged
				print("ADMM LOOP::convergence reached")
				if(sigma_threshold - maximum_frobenius_norm > 0.05*sigma_threshold) then--check if fake convergence
					print("ADMM LOOP::t'was fake convergence")
					--print("ADMM LOOP::scaling is: scaling= "..scaling)
					scaling = scaling*2.0;
					admm_steps=0
					--VecScaleAssign(u,1.0,u_converged);
					Jprime_ElemDisc:set_step_length(scaling);Jprime_DomainDisc:assemble_defect(SensitivityGF,u)
					SetZeroAwayFromSubset(SensitivityGF,"u1,u2,u3","obstacle_surface")
					print("ADMM LOOP::scaling is: scaling= "..scaling)
				else
					print("ADMM LOOP::convergence break case")
					break;
				end
				print("ADMM LOOP::scaling is: scaling= "..scaling)
			end	
			if admm_steps==admmSteps then
				print("ADMM LOOP::REACHED MAX ADMM STEPS");
				print("STEP "..step.." WILL BE REPEATED");
				p_solver_failure=true;admm_steps=0;
				break;	
			end																											
			admm_steps=admm_steps+1
	
		end--end of admm loop
		

		if bNewtonOutput then 
			gnuplot.write_data("__NewtonStats_step_"..step.."_.txt", {vNS_Step,vNS_NormSum,vNS_NormDeltaUp,
																												vNS_NormDeltaLambda,vNS_Lu2Norm}) 
			gnuplot.write_data("__NewtonIterations_step_"..step.."_.txt", {vNS_Step,vNS_RHS_Solver_Iterations,
								 			vNS_BvolSolver_Iterations,vNS_BxSolver_Iterations,vNS_BySolver_Iterations,vNS_LargeSolver_Iterations}) 
		end
		--check if solver failed, else check if a descent direction
		if p_solver_failure == true then
			steps_without_line_search = 0
			Lambda[1][1]=0.0
			Lambda[2][1]=0.0
			Lambda[3][1]=0.0
            Lambda[4][1]=0.0

			--step size control stuff
			scaling = 0.5*scaling
			--discard current iteration count
			sum_total_iterations = 0
			sum_largesolver_iterations = 0
			sum_bvolsolver_iterations = 0
			sum_bxsolver_iterations = 0
			sum_bysolver_iterations = 0
			sum_rhssolver_iterations = 0
			sum_newtonsteps = 0			
		else
			--Update geometry
			TransformDomainByDisplacement(u, "u1,u2,u3")

			--Solve primal problem on temporary vector and check if descent direction
			NavierStokes_Solver:init(navier_Op)----TODO:navierLine2
			if NavierStokes_Solver:prepare(v_temp) == false then print("NavierStokes::FAILED AFTER MESH DEFORMATION");exit();end
			if NavierStokes_Solver:apply(v_temp) == false then print("NavierStokes::FAILED AFTER MESH DEFORMATION");exit();end
			
			--Compute objective function on new geometry
			local drag_current = 0.5 * visc * Drag(u_zeros,v_temp,"v1,v2,v3","u1,u2,u3","outer",3)
			ShapeDerivative = VecProd(SensitivityGF,u)--here J_prime is additive and u_p is consistent
			vShapeDerivative[step]=ShapeDerivative

			--check if drag is descent direction
			print("NEW DRAG: "..drag_current)
			print("OLD DRAG: "..drag_old)
			print("NEW -OLD DRAG: "..drag_current - drag_old)
			if ((drag_current-drag_old) > 0.0 )then
				print("STEP: "..step.." NOT A DESCENT DIRECTION")
				VecScaleAssign(u_negative,-1.0,u);
				TransformDomainByDisplacement(u_negative, "u1,u2,u3")
				steps_without_line_search = 0
				sigma_threshold = 0.5*sigma_threshold
				Lambda[1][1]=0.0
				Lambda[2][1]=0.0
				Lambda[3][1]=0.0
                Lambda[4][1]=0.0
				--discard current iteration count
				sum_total_iterations = 0
				sum_largesolver_iterations = 0
				sum_bvolsolver_iterations = 0
				sum_bxsolver_iterations = 0
				sum_bysolver_iterations = 0
				sum_rhssolver_iterations = 0
				sum_newtonsteps = 0
			else
				--if steps_without_line_search > 10 then
					--step_control = math.min(2.0*step_control, step_control_init)
				--end

				print("STEP: "..step.." IS A DESCENT DIRECTION") 		
				VecScaleAssign(v,1.0,v_temp)--store temporary Navier-Stokes solutions			
				vDrag[step+1] = drag_current;--save new reference			
				gnuplot.write_data("__Drag.txt", {vStep,vDrag},false)			
				if bOutputFlows then
					print("Flows are saved to '" .. flowOutputFile .. "'...")
					vtkWriter:clear_selection()
					vtkWriter:select_nodal(flowNames, nodalFlows)
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print(flowOutputFile, v, step+1,step+1, false)
				end
				if bOutputPressure then
					print("Pressure field is saved")
					--PRESSURE
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("p", "pressure")
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print("pressureOutputFile", v,step,step, false)
				end
				--grid_positions:set(0.0);SaveNodalPositions2GridFunction(grid_positions,"u1,u2")
				--SaveToFile(grid_positions,"nodal_positions_"..step+1)
				if bDebugNodalPositions then
					print("Nodal positions saved to vtk")
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2,u3","u")
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print("grid_positions", grid_positions,step+1,step+1, false)
				end
				if bOutputMesh then
				print("u final is saved")
				vtkWriter:clear_selection()
				vtkWriter:select_nodal("u1,u2,u3","u")
				--vtkWriter:select_all(false)	-- write all variables
				vtkWriter:print("u", u,step+1,step+1, false)
				end

				vTotalLinearIterations[step] = sum_largesolver_iterations + sum_bvolsolver_iterations + sum_bxsolver_iterations + sum_bysolver_iterations + sum_bzsolver_iterations + sum_rhssolver_iterations
				vLargeSolverIterations[step] = sum_largesolver_iterations
				vBvolSolverIterations[step] = sum_bvolsolver_iterations
				vBxSolverIterations[step] = sum_bxsolver_iterations
				vBySolverIterations[step] = sum_bysolver_iterations
                vBzSolverIterations[step] = sum_bzsolver_iterations
				vRHSSolver[step] = sum_rhssolver_iterations
				vNonLinearIterationsPerP[step] = sum_newtonsteps
				gnuplot.write_data("__Iterations_per_step.txt", {vStep, vNonLinearIterationsPerP, vTotalLinearIterations, vRHSSolver, 
											vBvolSolverIterations, vBxSolverIterations, vBySolverIterations, vBzSolverIterations,vLargeSolverIterations})
				drag_old = drag_current
				break;
			end
		end--end if not solver failure
		--no break if solver failure, we just start over
	
	end--end of step size control loop
	step = step + 1

	u1_gradient_global=GlobalGridFunctionGradientData(u,"u1")
	u2_gradient_global=GlobalGridFunctionGradientData(u,"u2")
    u3_gradient_global=GlobalGridFunctionGradientData(u,"u3")
	u1_value_global=GlobalGridFunctionNumberData(u,"u1")
	u2_value_global=GlobalGridFunctionNumberData(u,"u2")
    u3_value_global=GlobalGridFunctionNumberData(u,"u3")

	--Reset the global data from Navier Stokes
	v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
	v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
    v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
	v1_value_global=GlobalGridFunctionNumberData(v,"v1")
	v2_value_global=GlobalGridFunctionNumberData(v,"v2")
    v3_value_global=GlobalGridFunctionNumberData(v,"v3")

	---Pressure
	p_value_global=GlobalGridFunctionNumberData(v,"p")

	--Reset the global data from Adjoint Navier Stokes
	q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
	q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
    q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
	q1_value_global=GlobalGridFunctionNumberData(q,"q1")
	q2_value_global=GlobalGridFunctionNumberData(q,"q2")
	q3_value_global=GlobalGridFunctionNumberData(q,"q3")
    ---Pressure
	h_value_global=GlobalGridFunctionNumberData(q,"h")

    lambda00_global=GlobalGridFunctionNumberData(lambda_piecewise,"l1")
    lambda01_global=GlobalGridFunctionNumberData(lambda_piecewise,"l2")
    lambda02_global=GlobalGridFunctionNumberData(lambda_piecewise,"l3")
    lambda10_global=GlobalGridFunctionNumberData(lambda_piecewise,"l4")
    lambda11_global=GlobalGridFunctionNumberData(lambda_piecewise,"l5")
    lambda12_global=GlobalGridFunctionNumberData(lambda_piecewise,"l6")
    lambda20_global=GlobalGridFunctionNumberData(lambda_piecewise,"l7")
    lambda21_global=GlobalGridFunctionNumberData(lambda_piecewise,"l8")
    lambda22_global=GlobalGridFunctionNumberData(lambda_piecewise,"l9")
    
    
    --projected value
    qproj00_global=GlobalGridFunctionNumberData(q_projected,"l1")
    qproj01_global=GlobalGridFunctionNumberData(q_projected,"l2")
    qproj02_global=GlobalGridFunctionNumberData(q_projected,"l3")
    qproj10_global=GlobalGridFunctionNumberData(q_projected,"l4")
    qproj11_global=GlobalGridFunctionNumberData(q_projected,"l5")
    qproj12_global=GlobalGridFunctionNumberData(q_projected,"l6")
    qproj20_global=GlobalGridFunctionNumberData(q_projected,"l7")
    qproj21_global=GlobalGridFunctionNumberData(q_projected,"l8")
    qproj22_global=GlobalGridFunctionNumberData(q_projected,"l9")

	--Set anew in AdjointNavierStokes
	AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
    AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
	AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)
    AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)

	DeformationEquationHessian_ElemDisc:set_deformation_d1(u1_value_global)
	DeformationEquationHessian_ElemDisc:set_deformation_d2(u2_value_global)
    DeformationEquationHessian_ElemDisc:set_deformation_d3(u3_value_global)
	DeformationEquationHessian_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	DeformationEquationHessian_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    DeformationEquationHessian_ElemDisc:set_deformation_vector_d3(u3_gradient_global)


	DeformationEquationRHS_ElemDisc:set_deformation_d1(u1_value_global)
	DeformationEquationRHS_ElemDisc:set_deformation_d2(u2_value_global)
    DeformationEquationRHS_ElemDisc:set_deformation_d3(u3_value_global)
	DeformationEquationRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	DeformationEquationRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    DeformationEquationRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	DeformationEquationRHS_ElemDisc:set_lambda00(lambda00_global)
	DeformationEquationRHS_ElemDisc:set_lambda01(lambda01_global)
	DeformationEquationRHS_ElemDisc:set_lambda02(lambda02_global)
	DeformationEquationRHS_ElemDisc:set_lambda10(lambda10_global)
	DeformationEquationRHS_ElemDisc:set_lambda11(lambda11_global)
	DeformationEquationRHS_ElemDisc:set_lambda12(lambda12_global)
	DeformationEquationRHS_ElemDisc:set_lambda20(lambda20_global)
	DeformationEquationRHS_ElemDisc:set_lambda21(lambda21_global)
	DeformationEquationRHS_ElemDisc:set_lambda22(lambda22_global)
	DeformationEquationRHS_ElemDisc:set_q00(qproj00_global)
	DeformationEquationRHS_ElemDisc:set_q01(qproj01_global)
	DeformationEquationRHS_ElemDisc:set_q02(qproj02_global)
	DeformationEquationRHS_ElemDisc:set_q10(qproj10_global)
	DeformationEquationRHS_ElemDisc:set_q11(qproj11_global)
	DeformationEquationRHS_ElemDisc:set_q12(qproj12_global)
	DeformationEquationRHS_ElemDisc:set_q20(qproj20_global)
	DeformationEquationRHS_ElemDisc:set_q21(qproj21_global)
	DeformationEquationRHS_ElemDisc:set_q22(qproj22_global)

	--Set Imports
	DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d1(u1_value_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d2(u2_value_global)
    DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_d3(u3_value_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    DeformationEquationLargeProblemRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda00(lambda00_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda01(lambda01_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda02(lambda02_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda10(lambda10_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda11(lambda11_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda12(lambda12_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda20(lambda20_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda21(lambda21_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_lambda22(lambda22_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q00(qproj00_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q01(qproj01_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q02(qproj02_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q10(qproj10_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q11(qproj11_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q12(qproj12_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q20(qproj20_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q21(qproj21_global)
	DeformationEquationLargeProblemRHS_ElemDisc:set_q22(qproj22_global)


	Jprime_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	Jprime_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
    Jprime_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	Jprime_ElemDisc:set_velocity_d1(v1_value_global)
	Jprime_ElemDisc:set_velocity_d2(v2_value_global)
    Jprime_ElemDisc:set_velocity_d3(v3_value_global)
	Jprime_ElemDisc:set_pressure(p_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
    Jprime_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
    Jprime_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
	Jprime_ElemDisc:set_adjoint_pressure(h_value_global)

	BVolume_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	BVolume_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    BVolume_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	BVolume_ElemDisc:set_deformation_d1(u1_value_global)
	BVolume_ElemDisc:set_deformation_d2(u2_value_global)
    BVolume_ElemDisc:set_deformation_d3(u3_value_global)

	XBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	XBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    XBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	XBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	XBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
    XBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

	YBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	YBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    YBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	YBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	YBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
    YBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

    ZBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	ZBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    ZBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	ZBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	ZBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
    ZBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

	MassModel_ElemDisc:set_deformation_d1(u1_value_global)
	MassModel_ElemDisc:set_deformation_d2(u2_value_global)
    MassModel_ElemDisc:set_deformation_d3(u3_value_global)
	MassModel_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	MassModel_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    MassModel_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	MassModel_ElemDisc:set_lambda00(lambda00_global)
	MassModel_ElemDisc:set_lambda01(lambda01_global)
	MassModel_ElemDisc:set_lambda02(lambda02_global)
	MassModel_ElemDisc:set_lambda10(lambda10_global)
	MassModel_ElemDisc:set_lambda11(lambda11_global)
	MassModel_ElemDisc:set_lambda12(lambda12_global)
	MassModel_ElemDisc:set_lambda20(lambda20_global)
	MassModel_ElemDisc:set_lambda21(lambda21_global)
	MassModel_ElemDisc:set_lambda22(lambda22_global)

	LambdaUpdate_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	LambdaUpdate_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
    LambdaUpdate_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	LambdaUpdate_ElemDisc:set_qproj00(qproj00_global)
	LambdaUpdate_ElemDisc:set_qproj01(qproj01_global)
	LambdaUpdate_ElemDisc:set_qproj02(qproj02_global)
	LambdaUpdate_ElemDisc:set_qproj10(qproj10_global)
	LambdaUpdate_ElemDisc:set_qproj11(qproj11_global)
	LambdaUpdate_ElemDisc:set_qproj12(qproj12_global)
	LambdaUpdate_ElemDisc:set_qproj20(qproj20_global)
	LambdaUpdate_ElemDisc:set_qproj21(qproj21_global)
	LambdaUpdate_ElemDisc:set_qproj22(qproj22_global)

end--end of optimization loop


