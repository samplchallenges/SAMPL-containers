# Uploading Your Challenge to [app.samplchallenges.org](https://app.samplchallenges.org/)

## Section 1: Pushing your container to a container repository
> Ahead of creating a submission to [app.samplchallenges.org](https://app.samplchallenges.org/), you will need to push your container image to any public container repository OR our provided private container repository. 
1. Use the `docker push <repository>/<username>/<container_name>:<tag>` command to push to your desired repository.
   * `repository` Examples: some popular repositories are listed below
     * Docker Hub: `docker.io`
     * GitHub Contaienr Repository: `ghcr.io`
   * `tag` Important information
     * When uploading your container to the repository for use on the [app.samplchallenges.org](https://app.samplchallenges.org/), it is *ESSENTIAL* that you use a unique tag every time you update and push your container.

## Section 2: Create a Submission on [app.samplchallenges.org](https://app.samplchallenges.org/)
1. Navigate to the [app.samplchallenges.org](https://app.samplchallenges.org/) website using a web browser of your choice.
2. Log in to [app.samplchallenges.org](https://app.samplchallenges.org/) using your credentials
3. Once at [app.samplchallenges.org](https://app.samplchallenges.org/), navigate to the list of Challenges using the "Challenges" tab at the top right or the "View Challenges" link. You will be redirected to the Challenges page. 
4. Click the name of the challenge you plan to participate in.
5. Click the "New Submission" button, which will open a New Submission Form. 
   * The Submission form is split into 4 sections: 
     * Container: Information about your container image
     * Special Arguments: Any extra command line keyword argument your container expects and the corresponding file. This section is intended for sensitive files and information, such as proprietary license files. 
     * Submission Details: Information about your submission such as methods, computing, etc.  
     * Submission Notes: Any annotations/notes-to-self about your submission. This section is for your benefit only. Challenge administrators will disregard any information in this section. This is the only section that can be modified once a challenge has ended.
7. In the Containers section, enter the information about your pushed container. This information should match the container information in your `docker push` command from Section 1
8. In the Special Arguments section, enter any extra keywork arguments and upload any corresponding files.
9. In the Submission Details section, enter the information about your submission, this section may be put in "draft mode" allowing you to come back and modify your submission details at a later time. 
10. Save your submission. You will be redirected to a page with a summary of your submission.

## Section 3: Run your submission.
1. 
