"""
pytest plugin file.

Adopts scenario-based test pattern from
https://docs.pytest.org/en/stable/example/parametrize.html#a-quick-port-of-testscenarios,
customising it to apply to test modules not test classes.
"""


def pytest_generate_tests(metafunc):
    """
    Parametrize tests using information about test scenarios.
    Test scenarios are assumed to be within a module-level
    ``SCENARIOS`` variable. This is assumed to be a dictionary
    of form::

        {
          <scenario-name>: {
            <fixture-name>-<fixture-value>
            [, <fixture-name>-<fixture-value>]*
          },
          ...
        }

    :param metafunc: pytest test function inspection object
    :type metafunc: _pytest.python.Metafunc
    """
    ids = []
    fix_values = []
    for scenario in metafunc.module.SCENARIOS:
        scenario_id, scenario_parameters = scenario
        ids.append(scenario_id)
        fix_params = scenario_parameters.items()
        fix_names = [fix_param[0] for fix_param in fix_params]
        fix_values.append(([fix_param[1] for fix_param in fix_params]))
    metafunc.parametrize(fix_names, fix_values, ids=ids, scope="module")
