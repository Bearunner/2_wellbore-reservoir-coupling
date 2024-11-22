//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"

class PorousFlowPeacemanBoreholeEnergyPro : public Kernel
{
public:
  static InputParameters validParams();

  PorousFlowPeacemanBoreholeEnergyPro(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const PorousFlowDictator & _dictator;
  const bool _var_is_porflow_var;
  const unsigned int _num_phases;

  const MaterialProperty<std::vector<Real>> & _pp;
  const MaterialProperty<std::vector<std::vector<Real>>> & _dpp_dvar;
  const MaterialProperty<Real> & _temperature;
  const MaterialProperty<std::vector<Real>> & _dtemperature_dvar;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _fluid_density_node;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables (at the qp)
  const MaterialProperty<std::vector<std::vector<Real>>> & _dfluid_density_node_dvar;
  /// PorousFlowDictator UserObject
  const MaterialProperty<std::vector<Real>> & _fluid_viscosity;

  /// Derivative of viscosity wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dfluid_viscosity_dvar;

  const MaterialProperty<std::vector<Real>> & _enthalpy;

  /// Derivative of the enthalpy wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _denthalpy_dvar;


const VariableValue & _p_well;
const VariableValue & _T_well;

const Real & _WI;
const Real & _dx;
const Real & _dy;
// const Real & _radius_ot;
// const Real & _lambda_oh;
const VariableValue & _rto;
const VariableValue & _Uto;
// const Real & _well_radius;
  /// Scale factor
const Real & _scale;
/// Optional function value
const Function & _function;
/// Optional Postprocessor value
const PostprocessorValue & _postprocessor;
const Real & _scale_T;
/// Optional function value
const Function & _function_T;
/// Optional Postprocessor value
const PostprocessorValue & _postprocessor_T;

Real computeQpJac(unsigned int pvar);
const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;

};
