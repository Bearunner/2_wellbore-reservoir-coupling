//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBoreholeMassS.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBoreholeMassS);

InputParameters
PorousFlowPeacemanBoreholeMassS::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<UserObjectName>(
       "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("p_well", "p_well");
  params.addRequiredCoupledVar("C_well", "C_well");
  params.addRequiredCoupledVar("mu_well", "mu_well");
  params.addParam<Real>("WI", "constant");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", 1, "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using PorousFlowPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

PorousFlowPeacemanBoreholeMassS::PorousFlowPeacemanBoreholeMassS(const InputParameters & parameters)
  : Kernel(parameters),
   _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
   _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
   _num_phases(_dictator.numPhases()),
   _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
   _dpp_dvar(
       getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
   _p_well(coupledValue("p_well")),
   _C_well(coupledValue("C_well")),
   _mu_well(coupledValue("mu_well")),
   _WI(this->template getParam<Real>("WI")),
   _scale(this->template getParam<Real>("value")),
   _function(getFunction("function")),
   _postprocessor(getPostprocessorValue("postprocessor"))

{
}

Real
PorousFlowPeacemanBoreholeMassS::computeQpResidual()

{
  Real r = 0.0;

  for (unsigned ph = 0; ph < _num_phases; ++ph)

    {
       r += _WI * _C_well[_qp] / _mu_well[_qp] * (_pp[_qp][ph] - _p_well[_qp]);
    }

  return r * _test[_i][_qp];
}


Real
PorousFlowPeacemanBoreholeMassS::computeQpJacobian()
{
  if (!_var_is_porflow_var)
  return 0.0;

  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}


Real
PorousFlowPeacemanBoreholeMassS::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}


Real
PorousFlowPeacemanBoreholeMassS::computeQpJac(unsigned int pvar)
{
Real j = 0.0;

for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
  j += _WI * _C_well[_qp] / _mu_well[_qp] * _dpp_dvar[_qp][ph][pvar] * _phi[_j][_qp];
  }

  return j * _test[_i][_qp];
}
