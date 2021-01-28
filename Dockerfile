FROM alpine:latest

ENV ATAT_FOLDER=/atatapp
ENV PATH "$PATH:/root/bin"

RUN apk update && apk add bash make g++ tcsh perl gnuplot

# RUN addgroup -S atatgroup && adduser -S atatuser -G atatgroup

COPY . /atatapp

#RUN mkdir -p atatapp && chown -R atatuser:atatgroup /atatapp && chmod -R 0755 /atatapp
#RUN mkdir -p atatapp
RUN chmod +x /atatapp/foolproof.sh

#USER atatuser
RUN mkdir -p /root/bin
WORKDIR /home/atatuser

RUN cd /atatapp && make && make install

CMD ["/bin/ash"]