"""
:py:mod:`riboviz.create_job_script` tests.
"""
import os
import tempfile
import pytest
import yaml
from riboviz import create_job_script
from riboviz import params


TEST_CONFIG = {
    params.JOB_NUM_CPUS: 1234,
    params.JOB_RUNTIME: "12:34:56",
    params.JOB_MEMORY: "2GB"
}
""" Test configuration. """
TEST_TEMPLATE = [
    create_job_script.TOKEN_TAG + key + create_job_script.TOKEN_TAG
    for key in TEST_CONFIG.keys()
]
""" Test template to be configured. """


@pytest.fixture(scope="function")
def tmp_config_file():
    """
    Create a temporary file with a ``yaml`` suffix and containing
    py:const:`TEST_CONFIG`.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    with open(tmp_file, 'w') as f:
        yaml.dump(TEST_CONFIG, f, default_flow_style=False)
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.fixture(scope="function")
def tmp_template_file():
    """
    Create a temporary file with a ``txt`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".txt")
    with open(tmp_file, 'w') as f:
        for line in TEST_TEMPLATE:
            f.write(line + '\n')
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.fixture(scope="function")
def tmp_output_file():
    """
    Create a temporary file with a ``txt`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".txt")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


def test_create_job_submission_script():
    """
    Test :py:func:`riboviz.create_job_script.create_job_submission_script`
    with a parameter for which there is a configuration value.
    """
    config = {params.JOB_NUM_CPUS: 123}
    template = ["prefix " +
                create_job_script.TOKEN_TAG +
                params.JOB_NUM_CPUS +
                create_job_script.TOKEN_TAG +
                " suffix"]
    custom = create_job_script.create_job_submission_script(config, template)
    assert custom == ["prefix 123 suffix"]


def test_create_job_submission_script_no_config():
    """
    Test :py:func:`riboviz.create_job_script.create_job_submission_script`
    with a parameter for which there is no configuration value.
    """
    config = {}
    template = ["prefix " +
                create_job_script.TOKEN_TAG +
                params.JOB_NUM_CPUS +
                create_job_script.TOKEN_TAG +
                " suffix"]
    custom = create_job_script.create_job_submission_script(config, template)
    assert custom == template


def test_create_job_submission_script_custom_token_tag():
    """
    Test :py:func:`riboviz.create_job_script.create_job_submission_script`
    with a parameter for which there is a configuration value and a
    custom token tag.
    """
    tag = "!!"
    config = {params.JOB_NUM_CPUS: 123}
    template = ["prefix " + tag + params.JOB_NUM_CPUS + tag + " suffix"]
    custom = create_job_script.create_job_submission_script(
        config, template, tag)
    assert custom == ["prefix 123 suffix"]


@pytest.mark.parametrize("value_flag",
                         [(True, "--" + params.VALIDATE_ONLY),
                          (False, "")])
def test_create_job_submission_script_validate_only(value_flag):
    """
    Test :py:func:`riboviz.create_job_script.create_job_submission_script`
    handling of :py:const:`riboviz.params.VALIDATE_ONLY` parameter.

    :param value_flag: Tuple with (value for
    :py:const:`riboviz.params.VALIDATE_ONLY`, expected flag in
    customised template)
    :type value_flag: tuple(bool, str or unicode)
    """
    value, flag = value_flag
    config = {params.VALIDATE_ONLY: value}
    template = ["prefix " +
                create_job_script.TOKEN_TAG +
                params.VALIDATE_ONLY +
                create_job_script.TOKEN_TAG +
                " suffix"]
    custom = create_job_script.create_job_submission_script(config, template)
    assert custom == ["prefix " + flag + " suffix"]


@pytest.mark.parametrize("value_prefix",
                         [("someone@ed.ac.uk", "#$ -M"),
                          (None, "# $ -M")])
def test_create_job_submission_script_job_email(value_prefix):
    """
    Test :py:func:`riboviz.create_job_script.create_job_submission_script`
    handling of :py:const:`riboviz.params.JOB_EMAIL` parameter in a
    line prefixed by ``#$ -M ...``.

    :param value_prefix: Tuple with (value for
    :py:const:`riboviz.params.JOB_EMAIL`, expected line prefix in
    customised template)
    :type value_prefix: tuple(bool, str or unicode)
    """
    value, prefix = value_prefix
    config = {params.JOB_EMAIL: value}
    template = ["#$ -M " +
                create_job_script.TOKEN_TAG +
                params.JOB_EMAIL +
                create_job_script.TOKEN_TAG]
    custom = create_job_script.create_job_submission_script(config, template)
    assert custom == [prefix + " " + str(value)]


@pytest.mark.parametrize("template_file", ["no-such-file.txt", os.getcwd()])
def test_create_job_script_config_file_error(tmp_config_file, template_file):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with missing or invalid template files.

    :param tmp_config_file: Test configuration file
    :type tmp_config_file: str or unicode
    :param template_file: Missing/invalid template file
    :type template_file: str or unicode
    """
    with pytest.raises(AssertionError):
        create_job_script.create_job_script(tmp_config_file, {}, template_file)


@pytest.mark.parametrize("config_file", ["no-such-file.txt", os.getcwd()])
def test_create_job_script_template_file_error(config_file, tmp_template_file):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with missing or invalid configuration files.

    :param config_file: Missing/invalid configuration file
    :type config_file: str or unicode
    :param tmp_template_file: Test template file
    :type tmp_template_file: str or unicode
    """
    with pytest.raises(AssertionError):
        create_job_script.create_job_script(config_file, {}, tmp_template_file)


def test_create_job_script_default(
        tmp_config_file, tmp_template_file, tmp_output_file):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with an empty configuration file and check that the template
    file has the default job configuration values applied.

    :param tmp_config_file: Test configuration file
    :type tmp_config_file: str or unicode
    :param tmp_template_file: Test template file
    :type tmp_template_file: str or unicode
    :param tmp_output_file: Test output file
    :type tmp_output_file: str or unicode
    """
    with open(tmp_config_file, 'w') as f:
        yaml.dump({}, f, default_flow_style=False)
    create_job_script.create_job_script(
        tmp_config_file, {}, tmp_template_file, tmp_output_file)
    with open(tmp_output_file) as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines]
    expected = [str(params.DEFAULT_JOB_CONFIG[key])
                for key in TEST_CONFIG.keys()]
    assert expected == lines


def test_create_job_script_config_file(
        tmp_config_file, tmp_template_file, tmp_output_file):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with a configuration file and check that the template
    file has the configuration values from the configuration
    file applied.

    :param tmp_config_file: Test configuration file
    :type tmp_config_file: str or unicode
    :param tmp_template_file: Test template file
    :type tmp_template_file: str or unicode
    :param tmp_output_file: Test output file
    :type tmp_output_file: str or unicode
    """
    create_job_script.create_job_script(
        tmp_config_file, {}, tmp_template_file, tmp_output_file)
    with open(tmp_output_file) as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines]
    expected = [str(TEST_CONFIG[key]) for key in TEST_CONFIG.keys()]
    assert expected == lines


def test_create_job_script_override(
        tmp_config_file, tmp_template_file, tmp_output_file):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with a configuration file and additional values and check that
    the template file has the additional values applied where
    they replace those in the configuration file.

    :param tmp_config_file: Test configuration file
    :type tmp_config_file: str or unicode
    :param tmp_template_file: Test template file
    :type tmp_template_file: str or unicode
    :param tmp_output_file: Test output file
    :type tmp_output_file: str or unicode
    """
    config = TEST_CONFIG.copy()
    config[params.JOB_NUM_CPUS] = 5678
    config[params.JOB_MEMORY] = "1TB"
    expected_config = config.copy()
    del config[params.JOB_RUNTIME]
    create_job_script.create_job_script(
        tmp_config_file, config, tmp_template_file, tmp_output_file)
    with open(tmp_output_file) as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines]
    expected = [str(expected_config[key]) for key in expected_config.keys()]
    assert expected == lines


def test_create_job_script_default_stdout(
        tmp_config_file, tmp_template_file, capsys):
    """
    Test :py:func:`riboviz.create_job_script.create_job_script`
    with an empty configuration file and check that the template
    file has the default job configuration values applied.

    :param tmp_config_file: Test configuration file
    :type tmp_config_file: str or unicode
    :param tmp_template_file: Test template file
    :type tmp_template_file: str or unicode
    """
    with open(tmp_config_file, 'w') as f:
        yaml.dump({}, f, default_flow_style=False)
    create_job_script.create_job_script(
        tmp_config_file, {}, tmp_template_file)
    captured = capsys.readouterr()
    lines = [line.strip() for line in captured.out.strip().split("\n")]
    expected = [str(params.DEFAULT_JOB_CONFIG[key])
                for key in TEST_CONFIG.keys()]
    assert expected == lines
