#
# www.htslib.org Deployment/Development Dockerfile
#
# https://github.com/samtools/www.htslib.org

FROM ubuntu:16.04
MAINTAINER "Joshua C. Randall" <jcrandall@alum.mit.edu>

# Install Prerequisites
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -q -y supervisor apache2 build-essential ruby-dev nodejs
RUN gem install jekyll jekyll-sitemap

# Add site source
ADD . /docker/www.htslib.org
WORKDIR /docker/www.htslib.org

# Build site with jekyll
RUN jekyll build --destination /var/www/html

# Add supervisord config
ADD _etc/docker-apache2-supervisord.conf /etc/supervisor/conf.d/apache2-supervisord.conf

# Start supervisord which will start apache2
ENTRYPOINT ["bash", "-c", "/usr/bin/supervisord -c /etc/supervisor/supervisord.conf ; $0 $*"]

# Default command watches the apache logs
CMD ["tail", "-f", "/var/log/apache2/*.log"]

