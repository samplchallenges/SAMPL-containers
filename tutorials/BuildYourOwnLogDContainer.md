# Building Your Own LogD Container
> This document details requirements and tips for writing a LogD container for SAMPL-challenges.

## Input Requirements
Every container must be able to handle the following input flags. These are the only flags your container will be expected to handle. We typically use [`click`](https://click.palletsprojects.com/en/8.0.x/) to handle command line argument parsing, but feel free to use your preferred parser.


## Output Requirements
**Printed Outputs**: Print the following to to `stdout`
The LAST line your container should output is below in the format key value where the key is logd and the value is the float LogD value. The key and value should be separated by a single space. You may print other outputs throughout your program, but this line must be the LAST line printed by your program.
```
logd {logd_float_value}
```

