# License and Code Privacy

We have included some options for participants who do not want their code and/or software licenses leaked or uploaded to a public repository. 

## Code Privacy
We realize that SAMPL participants may be concerned about sharing their code and other intellectual property. 

In order for the [app.samplchallenges.org](https://app.samplchallenges.org/) to have access to your containers upon submission we have provided two options:
* **Option 1:** Upload your container image YOUR OWN container repository ([ghcr.io](https://github.com/features/packages), [docker.io](https://www.docker.com/), etc.) and mark the image as PUBLIC. With this option, anyone can pull your container image and access any files you copied into the container image at build time.
* **Option 2:** Upload your container image to OUR PRIVATE container repository on Amazon Web Services. Our private container repository is a push only repository. Participants may only push NOT pull from this repo, meaning only [app.samplchallenges.org](https://app.samplchallenges.org/) website administrators from the Mobley Lab will have access to your container images and code. (Currently under construction)

To push to our private container repository on Amazon Web Services please do the following:
1. Add in steps for uploading to repository and discuss private key needed
   * Each time you push to the container repository, please use a UNIQUE tag for clarity. If you push to the container repository twice with the exact same tag, [app.samplchallenges.org](https://app.samplchallenges.org/) will not update your container to the newest version with that tag. You must use a new tag. 
      * For example, let's say I push a new container image as `ghcr.io/osatom/oedock:0.1`, create a submission with this image, and then realize I made a very small mistake and need to upload a modified version. When I push the modified container image, I would push with a new unique tag such as `0.1fix1`, written fully as `ghcr.io/osatom/oedock:0.1fix1`, rather than `ghcr.io/osatom/oedock:0.1`. 

## Proprietary Property and License Privacy
We realize that SAMPL participants will be concerned about leaking proprietary property and licenses of their sponsors or collaborators. 

To ensure only you and [app.samplchallenges.org](https://app.samplchallenges.org/) administrators have access to your license files, we have added a "Special Arguments" section to the submission page. Using the special arguments, you can add extra file keyword arguments to be called with your container at run time. Any keyword arguments you add to the "Special Arguments" section will be called with your container, meaning you should configure your container to accept the keyword argument and to expect its value to be a file. 
  * **Example**: In the case of a license, you can add a `license`, and upload your license file to be paired with they key. You also must ensure your container takes a `license` flag that expects a license file. The files you upload in the Special Arguments section are uploaded to a secure database and only accessible to you and website administrators; no other participants or website visitors will be able to access or use your license files. Your container will be called with `python ever_given/run.py <container_name> --license <license_file> --challenge_specific_args ...`
