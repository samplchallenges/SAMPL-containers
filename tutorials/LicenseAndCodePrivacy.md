# License and Code Privacy

We have included some options for participants who do not want their code and/or software licenses leaked or uploaded to a public repository. 

## Code Privacy
In order for the [app.samplchallenges.org](https://app.samplchallenges.org/) to have access to your containers upon submissionn there are two options:
1. your container can be uploaded to your container repository (ghcr.io, docker.io, etc.) and marked as public. With this option, 
2. your container can uploaded to our private container repository on Amazon Web Services. Our private container repository is a push only repository. Participants may only push NOT pull from it, making it so only app.samplchallenges.org website administrators from the Mobley Lab have access to your code.  

To push to our private container repository on Amazon Web Services please do the following:
Each time you push to the container repository, please use a UNIQUE tag for clarity. If you push to the container repository twice with the exact same tag, [app.samplchallenges.org](https://app.samplchallenges.org/) will not update your container to the newest version with that tag. You must use a new tag. 
* For example, let's say I push a new container image as `ghcr.io/osatom/oedock:0.1`, create a submission with this image, and then realize I made a very small mistake and need to upload a modified version. When I push the modified container image, I would push it as `ghcr.io/osatom/oedock:0.1fix1`, rather than `ghcr.io/osatom/oedock:0.1`. 
