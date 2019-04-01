FROM kernsuite/base:dev 

RUN docker-apt-install  casalite python-numpy python-casacore \ 
			meqtrees \
			python-pyxis \
			python-scatterbrane \			
			time\
			vim

#took out simms because installing it manually below

RUN pip install --upgrade pip

ENV LD_LIBRARY_CONFIG=/usr/local/lib

#ADD downloaded_files /downloaded_files
#RUN tar -xzf /downloaded_files/boost_1_64_0.tar.gz --directory / && rm /downloaded_files/boost_1_64_0.tar.gz
#RUN tar -xzf /downloaded_files/aatm-0.5.tar.gz --directory / && rm /downloaded_files/aatm-0.5.tar.gz
#RUN tar -xzf /downloaded_files/MeqSilhouette.tar.gz --directory / && rm /downloaded_files/MeqSilhouette.tar.gz

ADD downloaded_files/boost_1_64_0.tar.gz /
RUN cd /boost_1_64_0 && ./bootstrap.sh --with-libraries=program_options && ./b2 --prefix=/usr/local install

ADD downloaded_files/aatm-0.5.tar.gz /
RUN cd /aatm-0.5 && ./configure && make && make install

ADD downloaded_files/simms.tar.gz /
RUN cd /simms && python setup.py install

RUN ldconfig

RUN pip install termcolor 

ARG BUID=1001
RUN useradd -l -m -s /bin/bash -N -u $BUID mequser

ENV MEQTREES_CATTERY_PATH=/usr/lib/python2.7/dist-packages/Cattery/

ENV CODE_DIR=/vlbi-sim/
ADD /downloaded_files/vlbi-sim.tar.gz ${CODE_DIR}/..
RUN cd ${CODE_DIR}/framework/ && pip install -e .
RUN ln -sf $MEQTREES_CATTERY_PATH/Siamese/turbo-sim.py ${CODE_DIR}/framework/framework/turbo-sim.py 

RUN chown -R mequser ${CODE_DIR} 
USER mequser

# to complete initialization:
RUN echo -e "\nexit" | casa 
RUN python -c "import matplotlib.pyplot"

WORKDIR ${CODE_DIR}
CMD ["/bin/bash"] 
