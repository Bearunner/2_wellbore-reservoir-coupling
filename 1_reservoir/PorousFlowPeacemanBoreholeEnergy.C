//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBoreholeEnergy.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBoreholeEnergy);

InputParameters
PorousFlowPeacemanBoreholeEnergy::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<UserObjectName>(
       "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("p_well", "p_well");
  params.addRequiredCoupledVar("T_well", "T_well");
  params.addRequiredCoupledVar("rho_well", "rho_well");
  params.addRequiredCoupledVar("mu_well", "mu_well");
  params.addRequiredCoupledVar("h_well", "h_well");
  params.addParam<Real>("WI", "constant");
  params.addParam<Real>("dx", "dx");
  params.addParam<Real>("dy", "dy");
  // params.addParam<Real>("radius_ot", "outside radius of tubing");
  // params.addParam<Real>("lambda_oh", "thermal resistivity of openhole section");
  params.addRequiredCoupledVar("rto", "rto");
  params.addRequiredCoupledVar("Uto", "Uto");
  // params.addParam<Real>("well_radius", "constant");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", 1, "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.addParam<Real>("value_T", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function_T", 1, "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor_T", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  params.declareControllable("value_T");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using PorousFlowPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

PorousFlowPeacemanBoreholeEnergy::PorousFlowPeacemanBoreholeEnergy(const InputParameters & parameters)
  : Kernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _num_phases(_dictator.numPhases()),
    _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _dpp_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
    _temperature(getMaterialProperty<Real>("PorousFlow_temperature_qp")),
    _dtemperature_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_qp_dvar")),
    _p_well(coupledValue("p_well")),
    _T_well(coupledValue("T_well")),
    _rho_well(coupledValue("rho_well")),
    _mu_well(coupledValue("mu_well")),
    _h_well(coupledValue("h_well")),
    _WI(this->template getParam<Real>("WI")),
    _dx(this->template getParam<Real>("dx")),
    _dy(this->template getParam<Real>("dy")),
    // _radius_ot(this->template getParam<Real>("radius_ot")),
    // _lambda_oh(this->template getParam<Real>("lambda_oh")),
    _rto(coupledValue("rto")),
    _Uto(coupledValue("Uto")),
    // _well_radius(this->template getParam<Real>("well_radius")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _scale_T(this->template getParam<Real>("value")),
    _function_T(getFunction("function")),
    _postprocessor_T(getPostprocessorValue("postprocessor"))

{
}

Real
PorousFlowPeacemanBoreholeEnergy::computeQpResidual()

{
  Real r = 0.0;

  for (unsigned ph = 0; ph < _num_phases; ++ph)

  {
    r += _WI * (_h_well[_qp] + _p_well[_qp] / _rho_well[_qp]) * _rho_well[_qp] / _mu_well[_qp] * (_pp[_qp][ph] - _p_well[_qp]);

  //kinetic term
  // r += std::pow(_well_radius, 2)/8 * std::pow(_WI, 3) * _rho_well[_qp] / std::pow(_mu_well[_qp], 3) * std::pow((_pp[_qp][ph] - _p_well[_qp]), 3);

  // conduction term
   r += 2.0 * PI * _rto[_qp] * _Uto[_qp] * (_temperature[_qp] - _T_well[_qp]) / (_dx * _dy);
  }

  return r * _test[_i][_qp];
}


Real
PorousFlowPeacemanBoreholeEnergy::computeQpJacobian()
{
  if (!_var_is_porflow_var)
  return 0.0;

  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}


Real
PorousFlowPeacemanBoreholeEnergy::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}


Real
PorousFlowPeacemanBoreholeEnergy::computeQpJac(unsigned int pvar)
{
Real j = 0.0;

for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
  j += _WI * (_h_well[_qp] + _p_well[_qp] / _rho_well[_qp]) * _rho_well[_qp] / _mu_well[_qp] * _dpp_dvar[_qp][ph][pvar] * _phi[_j][_qp];
  j += 2.0 * PI * _rto[_qp] * _Uto[_qp] * _dtemperature_dvar[_qp][pvar] * _phi[_j][_qp] / (_dx * _dy);
  }

  return j * _test[_i][_qp];
}
