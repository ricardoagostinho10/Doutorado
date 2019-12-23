FROM python:onbuild
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt
EXPOSE 5000:5000
CMD python ./main.py
