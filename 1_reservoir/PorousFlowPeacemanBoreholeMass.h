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

class PorousFlowPeacemanBoreholeMass : public Kernel
{
public:

  static InputParameters validParams();

  PorousFlowPeacemanBoreholeMass(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // const unsigned int _fluid_component;
  const PorousFlowDictator & _dictator;
  const bool _var_is_porflow_var;
  const unsigned int _num_phases;

const MaterialProperty<std::vector<Real>> & _pp;
const MaterialProperty<std::vector<std::vector<Real>>> & _dpp_dvar;
const VariableValue & _p_well;
const VariableValue & _rho_well;
const VariableValue & _mu_well;
// MooseEnum _well_type;
const Real & _WI;
  /// Scale factor
const Real & _scale;
/// Optional function value
const Function & _function;
/// Optional Postprocessor value
const PostprocessorValue & _postprocessor;

//   /// Mass fraction of each component in each phase
 // const MaterialProperty<std::vector<std::vector<Real>>> & _mass_frac;
 // const MaterialProperty<std::vector<std::vector<std::vector<Real>>>> & _dmass_frac_dvar;

 Real computeQpJac(unsigned int pvar);
};
