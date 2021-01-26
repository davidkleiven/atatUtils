FROM alpine:latest

ENV ATAT_FOLDER=/atatapp
ENV PATH "$PATH:/home/atatuser/bin"

RUN apk update && apk add bash make g++ tcsh perl gnuplot

RUN addgroup -S atatgroup && adduser -S atatuser -G atatgroup

COPY . /atatapp

RUN mkdir -p atatapp && chown -R atatuser:atatgroup /atatapp && chmod -R 0755 /atatapp
RUN chmod +x /atatapp/foolproof.sh

USER atatuser

WORKDIR /home/atatuser
RUN mkdir -p /home/atatuser/bin
RUN cd /atatapp && make && make install

CMD ["/bin/ash"]