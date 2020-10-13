"""
Create a job submission script using values from a workflow
configuration file and a template job submission script.
"""
import os
import os.path
import yaml
from riboviz import params


def create_job_submission_script(config_file, template):
    """
    Create a job submission script using values from a workflow
    configuration and a template job submission script.

    :param config: Workflow configuration file
    :type config: str or unicode
    :param template: Job submission script template
    :type template: list(str or unicode)
    :return: Job submission script
    :rtype: list(str or unicode)
    """
    # TODO
    # Iterate through config_file
    # Search each linke for regexp token pattern e.g. %%PARAM%%
    # Replace marked up <PARAMETER> values with values
    # Look for tokens denoting additional content.
    print(config_file)
    print(template)
    return template


def create_job_script(config_file, input_file, output_file=None):
    """
    Create a job submission script using values from a workflow
    configuration file and a template job submission script.

    If ``output_file`` is ``None`` then the job submission script is
    printed to standard output.

    :param config_file: Workflow configuration file
    :type config_file: str or unicode
    :param input_file: Job submission script template
    :type input_file: str or unicode
    :param output_file: Job submission script
    :type output_file: str or unicode

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
    script = create_job_submission_script(config, template)
    if output_file is not None:
        with open(output_file, 'w') as f:
            f.writelines(["%s\n" % line  for line in script])
    else:
        for line in script:
            print(line)
