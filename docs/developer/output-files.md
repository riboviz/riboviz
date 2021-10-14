# Adding, renaming, and removing temporary or output files

* [Adding temporary or output files](#adding-temporary-or-output-files)
* [Renaming temporary or output files](#renaming-temporary-or-output-files)
* [Removing temporary or output files](#removing-temporary-or-output-files)

---

## Adding temporary or output files

If you add a workflow step that creates a new temporary or output file, then **respect the riboviz file naming conventions**, see [Input and output file names](style.md#input-and-output-file-names) and [Nextflow](style.md#nextflow) in [Style](./style.md).

You should also make the following updates.

### 1. Update Python code

If the file is produced by a riboviz Python script or a third-party tool then update `riboviz/workflow_files.py`, creating a Python constant for the new file.

If the file is produced by a riboviz R script then update `riboviz/workflow_r.py`, creating a Python constant for the new file.

Ensure the constant includes the file type (e.g. `TSV`, `PDF`, `H5`, `BAM` etc.) as its last component, delimited by an underscore.

For example:

```python
EXAMPLE_FILE_TSV = "example_file.tsv"
""" Example TSV file. """
```

If the file name will include a sample ID or other content which may change from run to run then define a Python format string instead, using `{}` for the part of the file name that may change. For example `riboviz/workflow_files.py` defines the following constants for sample-specific files:

```
ADAPTER_TRIM_FQ_FORMAT = "{}_trim.fq"

STATIC_HTML_FILE = "{}_output_report.html"
```

### 2. Update integration tests

Write a test function in `riboviz/test/integration/test_integration.py` to validate the temporary or output file by comparing it to a temporary or output file in the integration test data directory.

For information to help with writing an integration test, see [Writing an integration test](./testing.md#writing-an-integration-test.md) in [Developing and running tests](./testing.md).

### 3. Run integration tests

For the integration tests, tests will fail as the integration test data directory will not have the new file, so:

* Run the workflow.
* Copy the new file into the appropriate place in your integration test data directory.
* [Run vignette integration tests](./testing.md#run-vignette-integration-tests).

### 4. Update documentation

For temporary files:

* Update `docs/user/prep-riboviz-operation.md`, adding an entry for the new file to 'Temporary files'. Document any conditions which determine if the temporary file is produced e.g., dependencies on configuration parameters having certain values or the presence of certain input files.
* Update `docs/user/run-vignette.md`, adding the new file to 'Intermediate outputs in `vignette/tmp`'.

For output files:

* Update `docs/user/prep-riboviz-operation.md`, adding an entry for the new file to 'Output files'. Document any conditions which determine if the output file is produced e.g., dependencies on configuration parameters having certain values or the presence of certain input files.
 * Update `docs/user/riboviz-outputs.md`, adding a new section for the file with a description of the file, and add a link to this section from the links at the top of that page.
* Update `docs/user/run-vignette.md`, adding the new file to 'Outputs in `vignette/output`'.

For temporary and output files, update the workflow images, update the workflow images:

* Add the new file, with links from the process that creates it and to the process that consumes it, if applicable.
* Update the workflow SVG images, `docs/images/*.svg`.
* See [Updating workflow images](./documentation.md#updating-workflow-images) in [Writing and updating documentation](./documentation.md).

---

## Renaming temporary or output files

If renaming temporary or output files, then **respect the riboviz file naming conventions**, see [Input and output file names](style.md#input-and-output-file-names) and [Nextflow](style.md#nextflow) in [Style](./style.md).

You should also make the following updates.

### 1. Update Nextflow workflow

Search for all occurrences of the file name in `prep_riboviz.nf` and update these.

Rename any Nextflow variables that are derived from the file name, specifically channel names. For example, if renaming `old_example_file.tsv` to `new_example_file.tsv` and `prep_riboviz.nf` defines a channel `old_example_file_tsv`, then this channel would be renamed to `new_example_file_tsv`.

### 2. Update Python code

If the file is produced by a riboviz Python script or a third-party tool then update `riboviz/workflow_files.py`, renaming both the file and its associated Python constant.

If the file is produced by a riboviz R script then update `riboviz/workflow_r.py`, renaming both the file and its associated Python constant.

For example, if renaming `old_example_file.tsv` to `new_example_file.tsv` and `riboviz/workflow_files.py` defines:

```python
OLD_EXAMPLE_FILE_TSV = "old_example_file.tsv"
""" Example TSV file. """
```

then `riboviz/workflow_files.py` would be updated to:

```python
NEW_EXAMPLE_FILE_TSV = "new_example_file.tsv"
""" Example TSV file. """
```

Search for all references to the `riboviz/workflow_files.py` or `riboviz/workflow_r.py` constant as it was prior to renaming, across all the Python code, in `riboviz/` and its subdirectories, and update these references.

### 3. Run all tests

For the integration tests, tests will fail as the integration test data directory will not have the renamed file, so rename the file within your integration test data directory.

Run all tests:

* [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
* [Run R tests](./testing.md#run-r-tests).
* [Run vignette integration tests](./testing.md#run-vignette-integration-tests).

### 4. Update documentation

Search for all occurrences of the file name in the documentation, in `docs/` and its subdirectories, and update these.

Update the workflow images:

* Search for all occurrences of the file name in the workflow images, `docs/images/*.dot`, and update these.
* Update the workflow SVG images, `docs/images/*.svg`.
* See [Updating workflow images](./documentation.md#updating-workflow-images) in [Writing and updating documentation](./documentation.md).

---

## Removing temporary or output files

To remove a temporary or output file that is no longer created, make the following updates.

### 1. Update Nextflow workflow

Remove all occurrences of the file name and any associated channels from `prep_riboviz.nf`.

### 2. Update Python code

If the file is produced by a riboviz Python script or a third-party tool then update `riboviz/workflow_files.py`, removing both the file and its associated Python constant.

If the file is produced by a riboviz R script then update `riboviz/workflow_r.py`, removing both the file and its associated Python constant.

Search for all references to the `riboviz/workflow_files.py` or `riboviz/workflow_r.py` constant, across all the Python code, in `riboviz/` and its subdirectories, and remove these references.

### 3. Run all tests

Remove the file from your integration test data directory.

Run all tests:

* [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
* [Run R tests](./testing.md#run-r-tests).
* [Run vignette integration tests](./testing.md#run-vignette-integration-tests).

### 4. Update documentation

Search for all occurrences of the file name in the documentation, in `docs/` and its subdirectories, and remove these.

Update the workflow images:

* Search for all occurrences of the file name in the workflow images, `docs/images/*.dot`, and remove these.
* Update the workflow SVG images, `docs/images/*.svg`.
* See [Updating workflow images](./documentation.md#updating-workflow-images) in [Writing and updating documentation](./documentation.md).
