# Start from Ubuntu base image
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system tools and dependencies
RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    perl \
    wget \
    unzip \
    curl \
    bcftools \
    tabix \
    git \
    default-mysql-client \
    python3 \
    python3-pip \
    libdbi-perl \
    libdbd-mysql-perl \
    cpanminus \
    bwa \
    r-base \
    build-essential \
    libpq-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy Django app files
COPY . /app/

# Install Python dependencies
RUN pip3 install --upgrade pip && pip3 install --no-cache-dir -r requirements.txt

# --- Bioinformatics Tools ---

# Install Picard
RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar -O /opt/picard.jar

# Install GATK
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip -O /opt/gatk.zip \
    && unzip /opt/gatk.zip -d /opt/ && rm /opt/gatk.zip
ENV GATK_HOME=/opt/gatk-4.3.0.0
ENV PATH="${GATK_HOME}:${PATH}"

# Install SnpEff
RUN wget -q http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip && rm snpEff_latest_core.zip
ENV SNPEFF_HOME=/opt/snpEff
ENV PATH="${SNPEFF_HOME}:${PATH}"

# Install ANNOVAR
COPY bin/annovar.latest.tar /opt/
RUN tar -xvf /opt/annovar.latest.tar -C /opt/ && rm /opt/annovar.latest.tar
ENV ANNOVAR_HOME=/opt/annovar
ENV PATH="${ANNOVAR_HOME}:${PATH}"

# Set Picard jar path
ENV PICARD="/opt/picard.jar"

RUN mkdir -p /app/staticfiles
RUN chmod -R 755 /app/staticfiles

# Collect static files and run migrations
RUN python3 manage.py collectstatic --noinput
RUN python3 manage.py migrate --noinput
 
# Expose Django port
EXPOSE 8000

# Start server using gunicorn
CMD ["gunicorn", "visualization_workbench.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "3"]
