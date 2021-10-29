# Uploading Your Challenge to [app.samplchallenges.org](https://app.samplchallenges.org/)

## Section 1: Pushing your container to a container repository
> Ahead of creating a submission to [app.samplchallenges.org](https://app.samplchallenges.org/), you will need to push your container image to any public container repository OR our provided private container repository. 
1. Use the `docker push <repository>/<username>/<container_name>:<tag>` command to push to your desired repository.
   * `repository` Examples: some popular repositories are listed below
     * Docker Hub: `docker.io`
     * GitHub Contaienr Repository: `ghcr.io`
   * `tag` Important information
     * When uploading your container to the repository for use on the [app.samplchallenges.org](https://app.samplchallenges.org/), it is *ESSENTIAL* that you use a unique tag every time you update and push your container. For example, if I push a container called `docker.io/osatom/docking:0.1`, the next time I change and push that container I would push it as one of the following examples:
        * `docker.io/osatom/docking:0.1revised`
        * `docker.io/osatom/docking:0.1.1`
        * `docker.io/osatom/docking:0.2`

## Section 2: Create a Submission on [app.samplchallenges.org](https://app.samplchallenges.org/)
1. Navigate to the [app.samplchallenges.org](https://app.samplchallenges.org/) website using a web browser of your choice.
2. Log in to [app.samplchallenges.org](https://app.samplchallenges.org/) using your credentials
3. Once at [app.samplchallenges.org](https://app.samplchallenges.org/), navigate to the list of Challenges using the "Challenges" tab at the top right or the "View Challenges" link. You will be redirected to the Challenges page. 
4. Click the name of the challenge you plan to participate in.
5. Click the "New Submission" button, which will open a New Submission Form. 
   * The Submission form is split into 4 sections: 
     * **Container**: Information about your container image
     * **Special Arguments**: Any extra command line keyword argument your container expects and the corresponding file. This section is intended for sensitive files and information, such as proprietary license files. 
     * **Submission Details**: Information about your submission such as methods, computing, etc.  
     * **Submission Notes**: Any annotations/notes-to-self about your submission. This section is for your benefit only. Challenge administrators will disregard any information in this section. This is the only section that can be modified once a challenge has ended.
7. In the Containers section, enter the information about your pushed container. This information should match the container information in your `docker push` command from Section 1
   * If I used the command `docker push -t ghcr.io/osatom/docking:0.1revised`to push my container, the Container section fields of my Submission would look like the following:
      * **Name**: Docking Container - 0.1revised (Please note the container name can be anything you desire)
      * **Registry**: ghcr.io
      * **Label**: osatom/docking
      * **Tag**: 0.1revised
9. In the Special Arguments section, enter any extra keyword arguments and upload any corresponding files that are necessary to run your container
10. In the Submission Details section, enter the information about your submission, this section may be put in "draft mode" allowing incomplete fields to be saved and finished at a later time. If you are unsure of your submission details, please un-check the "Ranked" box. 
12. Save your submission. You will be redirected to a page with a summary of your submission.

## Section 3: Run your submission.
1. If you just finished Section 2, skip to step 5.
2. If you are returning to the website, ensure you are logged in, otherwise your submissions will not be visible.
3. Click on the "Submissions" tab at the top of the webpage.
4. Click on the name of the submission you want to run.
5. Click the blue "Run Submission" button. 
6. View the `stdout` or `stderr` output from your container by clicking on any bulleted link in the "Public Run" section to view the output.
7. Wait for your submission to finish running, this may take a while. Please note, if your container fails on the Public Dataset, the container will not be run with the Private Dataset.

## Section 4: Finalizing your submission
1. When you are happy with your submission, in the Submission Details section you will need to un-check "Draft mode" on your submission and check the "Ranked" box.

## Section 5: After the Challenge has ended
1. After the challenge has ended you may log in, and return to your submission details page.
2. If you'd like to leave notes related to your submission or its performance, click "Update" and type your notes in the "Submission Notes" section. Only this section will be modifiable post-challenge.



