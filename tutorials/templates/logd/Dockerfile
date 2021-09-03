  
FROM 
# inherit from base container FROM <base container>

WORKDIR /opt/app/

ENV PATH="/root/.local/bin:$PATH"

COPY setup.py main_template.py ./
# copy in files you need COPY <files> ./

RUN pip install .

RUN ls -l /opt/app

ENTRYPOINT ["run-logd"]
# add in entry_point from ENTRYPOINT ["<entrypoint_from_setup.py>"]
