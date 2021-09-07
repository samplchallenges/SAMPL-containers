FROM adv-tutorial-base

WORKDIR /opt/app/
# set the work directory

COPY setup.py autodock.py main.py ./

RUN pip install .
# install setup.py

ENV PATH="/root/.local/bin:$PATH"

RUN ls -l /opt/app

ENTRYPOINT ["run-autodock"]
