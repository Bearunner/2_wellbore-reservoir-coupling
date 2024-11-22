//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBoreholeMassPro.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBoreholeMassPro);

InputParameters
PorousFlowPeacemanBoreholeMassPro::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<unsigned int>(
    "fluid_component", 0, "The index corresponding to the component for this kernel");
  params.addRequiredParam<UserObjectName>(
       "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("p_well", "p_well");
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

PorousFlowPeacemanBoreholeMassPro::PorousFlowPeacemanBoreholeMassPro(const InputParameters & parameters)
  : Kernel(parameters),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _num_phases(_dictator.numPhases()),
    _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _dpp_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
    _fluid_density_node(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")),
    _dfluid_density_node_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_density_nodal_dvar")),
    _fluid_viscosity(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")),
    _dfluid_viscosity_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_nodal_dvar")),
    _p_well(coupledValue("p_well")),
    _WI(this->template getParam<Real>("WI")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _mass_frac(getMaterialProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal")),
    _dmass_frac_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real>>>>(
     "dPorousFlow_mass_frac_nodal_dvar"))

{
  if (_fluid_component >= _dictator.numComponents())
  paramError(
      "fluid_component",
      "The Dictator proclaims that the maximum fluid component index in this simulation is ",
      _dictator.numComponents() - 1,
      " whereas you have used ",
      _fluid_component,
      ". Remember that indexing starts at 0. The Dictator does not take such mistakes lightly.");
}

Real
PorousFlowPeacemanBoreholeMassPro::computeQpResidual()
{
    Real r = 0.0;

    for (unsigned ph = 0; ph < _num_phases; ++ph)
    {
     r += _WI * _fluid_density_node[_i][ph] / _fluid_viscosity[_i][ph] * (_pp[_qp][ph] - _p_well[_qp]) * _mass_frac[_i][ph][_fluid_component];
    }

  return r * _test[_i][_qp];
}


Real
PorousFlowPeacemanBoreholeMassPro::computeQpJacobian()
{
  if (!_var_is_porflow_var)
  return 0.0;

  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}


Real
PorousFlowPeacemanBoreholeMassPro::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}


Real
PorousFlowPeacemanBoreholeMassPro::computeQpJac(unsigned int pvar)
{
Real j = 0.0;

for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
  j += _WI *_fluid_density_node[_i][ph] / _fluid_viscosity[_i][ph] * _dpp_dvar[_qp][ph][pvar] * _phi[_j][_qp] * _mass_frac[_i][ph][_fluid_component];
  }

if (_i != _j)
  return j * _test[_i][_qp];

for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
  j += _WI * (_pp[_qp][ph] - _p_well[_qp]) * _mass_frac[_i][ph][_fluid_component] * (_dfluid_density_node_dvar[_i][ph][pvar] / _fluid_viscosity[_i][ph] - _dfluid_viscosity_dvar[_i][ph][pvar] * _fluid_density_node[_i][ph] / std::pow(_fluid_viscosity[_i][ph], 2));
  j += _WI * (_pp[_qp][ph] - _p_well[_qp]) * _dmass_frac_dvar[_i][ph][_fluid_component][pvar] * _fluid_density_node[_i][ph] / _fluid_viscosity[_i][ph];
  }

  return j * _test[_i][_qp];
}
