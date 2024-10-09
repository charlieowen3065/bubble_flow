#include "BubbleFlowApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
BubbleFlowApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

BubbleFlowApp::BubbleFlowApp(InputParameters parameters) : MooseApp(parameters)
{
  BubbleFlowApp::registerAll(_factory, _action_factory, _syntax);
}

BubbleFlowApp::~BubbleFlowApp() {}

void
BubbleFlowApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<BubbleFlowApp>(f, af, s);
  Registry::registerObjectsTo(f, {"BubbleFlowApp"});
  Registry::registerActionsTo(af, {"BubbleFlowApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
BubbleFlowApp::registerApps()
{
  registerApp(BubbleFlowApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
BubbleFlowApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  BubbleFlowApp::registerAll(f, af, s);
}
extern "C" void
BubbleFlowApp__registerApps()
{
  BubbleFlowApp::registerApps();
}
