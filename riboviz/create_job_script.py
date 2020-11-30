"""
Create a job submission script using values from a workflow
configuration file and a template job submission script.
"""
import os
import os.path
import re
import yaml
from riboviz import params


TOKEN_TAG = "%%"
"""
Default tag marking up tokens in template job submission script to be
replaced by configuration values.
"""
JOB_CONFIG = {
    params.JOB_NAME: "riboviz",
    params.JOB_RUNTIME: "48:00:00",
    params.JOB_MEMORY: "8GB",
    params.JOB_NUM_CPUS: 16,
    params.JOB_EMAIL_EVENTS: "beas",
    params.JOB_EMAIL: None,
    params.NEXTFLOW_WORK_DIR: "work",
    params.CREATE_NEXTFLOW_REPORT: False,
    params.NEXTFLOW_REPORT_FILE: "report.html",
    params.VALIDATE_ONLY: False
}
""" Default values for job configuration parameters. """


def create_job_submission_script(config,
                                 template,
                                 token_tag=TOKEN_TAG):
    """
    Create a job submission script using values from a workflow
    configuration and a template job submission script.

    :param config: Workflow configuration
    :type config: str or unicode
    :param template: Job submission script template
    :type template: list(str or unicode)
    :param token_tag: Tag marking up tokens in template job submission
    script to be replaced by configuration values.
    :type token_tag: str or unicode
    :return: Job submission script
    :rtype: list(str or unicode)
    """
    token_regex = "{tag}.*?{tag}".format(tag=token_tag)
    script = []
    for line in template:
        nu_line = line
        matches = re.finditer(token_regex, line)
        for match in matches:
            tag = match.group(0).strip(token_tag)
            if tag in config:
                nu_line = nu_line.replace(match.group(0), str(config[tag]))
        script.append(nu_line)
    return script


def create_job_script(config_file,
                      override_config,
                      input_file,
                      output_file=None,
                      token_tag=TOKEN_TAG):
    """
    Create a job submission script using values from a workflow
    configuration file and a template job submission script.

    If ``output_file`` is ``None`` then the job submission script is
    printed to standard output.

    :param config_file: Workflow configuration file
    :type config_file: str or unicode
    :param override_config: Extra configuration, overriding that in
    ``config_file`` (if there are any duplicates between these)
    :type override_config: dict(str or unicode => str or unicode)
    :param input_file: Job submission script template
    :type input_file: str or unicode
    :param output_file: Job submission script
    :type output_file: str or unicode
    :param token_tag: Tag marking up tokens in template job submission
    script to be replaced by configuration values.
    :type token_tag: str or unicode
    :raises AssertionError: If ``config_file`` or
    ``input_file`` does not exist or is not a file
    """
    for file_name in [input_file, config_file]:
        assert os.path.exists(file_name) and os.path.isfile(file_name),\
            "{} does not exist or is not a file".format(file_name)
    with open(input_file, 'r') as f:
        template = [line.rstrip() for line in f]
    with open(config_file, 'r') as f:
        file_config = yaml.load(f, yaml.SafeLoader)

    config = JOB_CONFIG.copy()
    config.update(file_config)
    config.update(override_config)
    script = create_job_submission_script(config, template, token_tag)
    if output_file is not None:
        with open(output_file, 'w') as f:
            f.writelines(["%s\n" % line for line in script])
    else:
        for line in script:
            print(line)
