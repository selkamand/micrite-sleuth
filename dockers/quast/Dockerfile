FROM python:3.11.15-alpine3.23

# Install basic system packages
RUN apk add bash

# Install software we're after 
RUN pip install "quast==5.2.0"

# Remove unnecessary packages to avoid vulnerabilities
RUN pip uninstall -y wheel setuptools
