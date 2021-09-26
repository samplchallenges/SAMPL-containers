# License and Code Privacy

We have included some options for participants who do not want their code and/or software licenses leaked or uploaded to a public repository. 

## Code Privacy
In order for the [app.samplchallenges.org](https://app.samplchallenges.org/) to have access to your containers upon submission, you containers must be uploaded to your container repository (ghcr.io, docker.io, etc.) and marked as public or uploaded to our private container repository on Amazon Web Services. Our private container repository is a push only repository. Participants may only push NOT pull from it, making it so only app.samplchallenges.org website administrators from the Mobley Lab have access to your code.  

To push to the private container repository please do the following:

Each time you push to the containe repository, please use a UNIQUE tag for clarity. If you push to the container repository twice with the exact same tag, [app.samplchallenges.org](https://app.samplchallenges.org/) will not update your container to the newest version with that tag. You must use a new tag. 
* For example, if I push a new container as `ghcr.io/osatom/oedock:0.1` and realize I made a small mistake and need to upload a modified version, I would push the new container as `ghcr.io/osatom/oedock:0.1fix1` rather than `ghcr.io/osatom/oedock:0.1` again. 
