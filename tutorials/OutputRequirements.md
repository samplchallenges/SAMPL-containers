# Output Requirements

## File Required Outputs: 
If a required output is a file, save the file to to directory passed to your container using the `--output-dir`. Your file name can be anything as long as you save the file to specified `--output-dir` and print out the full file path as a printed output with the required key.
* Example: `path_to_output_file = {output_dir}/{output_file_name}`

## Printed Required Outputs (`stdout`): 
Whether the required output is a numerical value or a file, you will ALWAYS need to print at least one output from your container in the format `key value` where the `key` is a string value specific to the challenge. The LAST lines your container should output to `stdout` are your printed required outputs. The key and value should be separated by a single space. You may print other outputs throughout your program, but the required printed outputs must be the LAST lines printed to `stdout` by your program.
* The container would output something similar to the following
    * Numerical Example: `key numerical_value`
    * File Example: `key path_to_output_file` 
         * Please see previous section "File Required Ouputs" for more details about `path_to_output_file`

## Intentional No Prediction Output: 
* Protocols for intentional no prediction for an input will be determined for each specific challenge. Please see the specific challenge for more details.

## Program Logs
* Any output to `stdout` or `stderr` will be logged with timestamps associated with each output. These logs will be made accessible to you.
* Please print general logging info to `stdout` and error messages to `stderr` as is convention.
* As stated in "Printed Required Outputs", the last two lines of `stdout` output must be your required `key value` pairs. 
