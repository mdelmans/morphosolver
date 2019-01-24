FROM quay.io/fenicsproject/stable:2018.1.0.r3

RUN mkdir -p /morphosolverhome
RUN chown fenics /morphosolverhome
WORKDIR /morphosolverhome

COPY morphosolver /

RUN git clone https://github.com/mikaem/fenicstools.git --branch 2018.1 /morphosolverhome/fenicstools

RUN sudo pip3 install cppimport==18.11.8

RUN sudo pip3 install /morphosolverhome/fenicstools/

RUN sudo python3 /morphosolverhome/fenicstools/tests/test_Interpolation.py

RUN sudo rm -rf /morphosolverhome/fenicstools

RUN sudo pip3 install git+https://github.com/mdelmans/morphosolver.git