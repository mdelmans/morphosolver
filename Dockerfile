FROM quay.io/fenicsproject/stable:current

RUN mkdir -p /morphosolverhome
RUN chown fenics /morphosolverhome
WORKDIR /morphosolverhome

COPY morphosolver /

RUN git clone https://github.com/mikaem/fenicstools.git /morphosolverhome/fenicstools

RUN sudo pip3 install cppimport

RUN sudo pip3 install /morphosolverhome/fenicstools/

RUN sudo python3 /morphosolverhome/fenicstools/tests/test_Interpolation.py

RUN sudo rm -rf /morphosolverhome/fenicstools

RUN sudo pip3 install git+https://github.com/mdelmans/morphosolver.git