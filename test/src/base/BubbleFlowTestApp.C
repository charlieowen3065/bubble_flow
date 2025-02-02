//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "BubbleFlowTestApp.h"
#include "BubbleFlowApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
BubbleFlowTestApp::validParams()
{
  InputParameters params = BubbleFlowApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

BubbleFlowTestApp::BubbleFlowTestApp(InputParameters parameters) : MooseApp(parameters)
{
  BubbleFlowTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

BubbleFlowTestApp::~BubbleFlowTestApp() {}

void
BubbleFlowTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  BubbleFlowApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"BubbleFlowTestApp"});
    Registry::registerActionsTo(af, {"BubbleFlowTestApp"});
  }
}

void
BubbleFlowTestApp::registerApps()
{
  registerApp(BubbleFlowApp);
  registerApp(BubbleFlowTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
BubbleFlowTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  BubbleFlowTestApp::registerAll(f, af, s);
}
extern "C" void
BubbleFlowTestApp__registerApps()
{
  BubbleFlowTestApp::registerApps();
}
