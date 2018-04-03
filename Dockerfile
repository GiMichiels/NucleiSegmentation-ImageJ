FROM neubiaswg5/fiji-base:latest

RUN cd plugins && wget -O imagescience.jar https://imagescience.org/meijering/software/download/imagescience.jar
RUN cd plugins && wget -O FeatureJ_.jar https://imagescience.org/meijering/software/download/FeatureJ_.jar

ADD macro.ijm /fiji/macros/macro.ijm                                           
ADD run.sh /fiji/run.sh
ADD wrapper.py /
RUN cd /fiji && chmod a+x run.sh

FROM gmichiels/python-client-base

ENTRYPOINT ["python", "/wrapper.py"]

