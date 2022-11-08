# Writing Fail-safe Containers
Due to the automated nature of the containerized challenges, containers will need to be as fail-safe as possible. In other words, your program should fail gracefully on input(s), and still output a default result. 

## Running processes:
To run processes or command line programs from your container, we recommend using the [`subprocess`](https://docs.python.org/3/library/subprocess.html) python library. Using the subprocess library allows you to:
1. Run processes and command line programs
   * `subprocess.run()`
2. Set a timeout for processes so your container will not run forever if it gets stuck on an input: 
   * `subprocess.run(timeout=[seconds])` 
3. Allows you to check the returncode of your finished process. If you have multiple steps in your workflow that require processes, this will allow you to check for errors earlier in your workflow
   * `subprocess.CompletedProcess.check_returncode()`


**A note about `subprocess.run(shell=True)`:**
* Oftentimes, the `shell=True` option for `subprocess.run()` is overused. 
* The shell used when `shell=True` is determined by the `SHELL` environment variable, which means you are potentially invoking a *mystery* program.
* Using `shell=True` launches a full shell (which is often unnecessary overhead) and adds more complexity to your program, meaning more places your program can fail 
* To avoid using `shell=True` you will need to tokenize your command
   * Instead of this: `subprocess.run("echo hello", shell=True)`
   * Use this: `subprocess.run(["echo", "hello"])`
* An informative discussion of how `shell=True` works can be found [here](https://discuss.dizzycoding.com/actual-meaning-of-shelltrue-in-subprocess/)

## Error catching:
If you know / suspect your program will raise errors, it is best to use [`try` and `except` and `finally`](https://docs.python.org/3/tutorial/errors.html#handling-exceptions) blocks to allow your program to fail on inputs gracefully.
* For example, if I used `subprocess.run(['run','long','program'], timeout=60*10)` I would expect a `subprocess.TimeoutExpired` Exception if my program surpassed 10 minutes (or 60 * 10 seconds)
   * To handle this I would write the following code:
      ```
      try:
          subprocess.run(['run','long','program'], timeout=60*10)
      except subprocess.TimeoutExpired:
          print('TimeoutExpired: max program run time surpassed')
          # output null or default result
      ```
* Because we would like our containers to fail gracefully and "catch all" failures I would do something like the following:
      ```
      # create a default result to output in case of failure
      result = -1 
      
      try:
          # try the following code:
          # program code goes here
          # use your program code to determine your result 
          result = 10000
      except Exception as e:
          print(f"ERROR: {e}")
      finally:
          # export out the result whether default or calculated
          # whether as a print statement (for un-batched) or as a text file (for batched)
      ```
      * A more in depth example can be found [here](https://github.com/samplchallenges/SAMPL-containers/blob/main/virtual_screening/examples/adv-screen-docker/main.py)

* For a nice tutorial on `try` and `except` blocks, see [geeksforgeeks](https://www.geeksforgeeks.org/python-try-except/)
