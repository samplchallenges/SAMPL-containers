# Uploading Your Container to SAMPL Private AWS Elastic Container Registry (ECR)
In this document, we detail how to upload your container to our private push-only AWS container registry. If you push your container to docker.io or ghcr.io, our infrastructure will only be able to access your container if it is marked as public. However, marking your container as public allows anyone to potentially access your code and/or softwares installed in your container. By pushing your container to our private push-only registry, your container is private to everyone except SAMPL challenge administrators and our infrastructure. This provides much more security for your code and/or softwares installed in your container.

## Set Up AWS Command Line Tool
1. If you haven't already, please contact the SAMPL Challenges administrators for the `AWS Access Key ID` and `AWS Secret Access Key`
2. Download the aws command line tool from the following link https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
3. Please choose from one of the two options to set up your credentials

### Option 1: Setup Credentials
1. From the command line, run the command: aws config
2. When prompted by the aws config command enter `AWS Access Key ID` and `AWS Secret Access Key` from the SAMPL Challenges administrators. You may leave `Default output format` blank by hitting enter.
    ```
    AWS Access Key ID [None]: <from administrators>
    AWS Secret Access Key [None]: <from administrators>
    us-east-2 [None]: us-east-2
    Default output format [None]: 
    ```
3. These are now your default AWS credentials, you will not need to specify a `profile` when you run the `aws login` command

### Option 2: Setup Credentials
1. Create a directory in your home directory called `.aws` 
    * `mkdir ~/.aws`
2. Change directories into `~/.aws`
3. Open a file for editing called `config`
4. In the config file, add the following:
  ```
  [sampl_push]
  region = us-east-2
  output = None
  ```
5. Save `config` file
6. Open a file for editing called `credentials`
7. In the `credentials` file, add the following, replacing `<from administrators>` with the `AWS Access Key ID` and `AWS Secret Access Key` from the SAMPL Challenges administrators. 
  ```
  [sampl_push]
  aws_access_key_id = <from administrators>
  aws_secret_access_key = <from administrators>
  ```
8. Save `credentials` file
9. These are now your AWS credentials under the `profile` called `sampl_push`, you will need to specify the `profile` when you run the `aws login` command


## Prepare Your Container for Upload
These instructions assume that you have already built your container using `docker build`.
1. Use the `docker images` command to get either the `IMAGE_ID` or the `image_name:tag` of the image you would like to upload
2. Determine the tag for your image. The TAG must use the pattern `[challenge_tag]-[first_initial_last_name]-[container_name_and_version]`. The `challenge_tag` will be specified on an on challenge basis. Please check the specific challenge page for instructions. If you re-upload a newer version of your container, please increment the version number as you cannot overwrite images. 
3. Re-tag the image you would like to upload using the following command where your IMAGE_ID is from step 1 and the TAG uses the pattern from step 2. 
  * `docker tag [IMAGE_ID] 103125031445.dkr.ecr.us-east-2.amazonaws.com/sampl-league:[TAG]`
6. Login to the aws repository 
  * If you used [Option 1: Setup Credentials]():
    * `aws ecr get-login-password --region us-east-2 | docker login --username AWS --password-stdin 103125031445.dkr.ecr.us-east-2.amazonaws.com`
  * If you used [Option 2: Setup Credentials]():
    * `aws ecr get-login-password --region us-east-2 --profile sampl_push | docker login --username AWS --password-stdin 103125031445.dkr.ecr.us-east-2.amazonaws.com`
7. Push your container to the AWS push only repository
  * `docker push 103125031445.dkr.ecr.us-east-2.amazonaws.com/sampl-league:[TAG]`
8. If you are concerned you can try to pull your image to ensure others cannot download it
  * `docker pull 103125031445.dkr.ecr.us-east-2.amazonaws.com/sampl-league:[TAG]`

## Create a Submission on app.samplchallenges.org
1. Select the challenge you are participating in.
2. Click the "New Submission" button or if you have already started a submission, navigate to your submission and click the "Update" button.
3. Navigate to the Container header
4. Fill out the Registry, Label, and Tag as follows
    ```
    Registry: 103125031445.dkr.ecr.us-east-2.amazonaws.com
    Label: sampl-league
    Tag: [TAG_from_step_2]
    ```
