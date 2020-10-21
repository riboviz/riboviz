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
                nu_line = nu_line.replace(match.group(0), config[tag])
        script.append(nu_line)
    return script


def create_job_script(config_file,
                      extra_config,
                      input_file,
                      output_file=None,
                      token_tag=TOKEN_TAG):
    """
    Create a job submission script using values from a workflow
    configuration file and a template job submission script.

    If ``output_file`` is ``None`` then the job submission script is
    printed to standard output.

    Two values are added to the configuration:

    * ``config_file``: base-name (local name) of ``config_file`.
    * ``config_path``: absolute path to ``config_file``, including the
      base-name.

    :param config_file: Workflow configuration file
    :type config_file: str or unicode
    :param extra_config: Extra configuration, overriding that in
    ``config_file`` (if there are any duplicates between these)
    :type extra_config: dict(str or unicode => str or unicode)
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
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    with open(input_file, 'r') as f:
        template = [line.rstrip() for line in f]

    config["config_file"] = os.path.basename(config_file)
    config["config_path"] = os.path.abspath(config_file)
    # Configuration values in extra_config override those in config.
    config.update(extra_config)
    script = create_job_submission_script(config,
                                          template,
                                          token_tag)
    if output_file is not None:
        with open(output_file, 'w') as f:
            f.writelines(["%s\n" % line for line in script])
    else:
        for line in script:
            print(line)
