# -----------------------------
# Base image
# -----------------------------
    FROM ubuntu:22.04

    # -----------------------------
    # Environment setup
    # -----------------------------
    ENV DEBIAN_FRONTEND=noninteractive
    ENV PYTHONDONTWRITEBYTECODE=1
    ENV PYTHONUNBUFFERED=1
    
    # Set working directory
    WORKDIR /app
    
    # -----------------------------
    # Install system packages
    # -----------------------------
    RUN apt-get update && apt-get install -y \
        build-essential \
        openjdk-11-jre-headless \
        perl \
        wget \
        unzip \
        curl \
        git \
        bcftools \
        tabix \
        default-mysql-client \
        libdbi-perl \
        libdbd-mysql-perl \
        cpanminus \
        bwa \
        r-base \
        libpq-dev \
        python3.11 \
        python3.11-venv \
        python3.11-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*
    
    # -----------------------------
    # Python environment
    # -----------------------------
    RUN python3.11 -m venv /opt/venv
    ENV PATH="/opt/venv/bin:$PATH"
    
    # Upgrade pip and install requirements
    COPY requirements.txt /app/
    RUN pip install --upgrade pip \
        && pip install --no-cache-dir -r requirements.txt
    
    # -----------------------------
    # Copy application code
    # -----------------------------
    COPY . /app/
    
    # -----------------------------
    # Bioinformatics tools
    # -----------------------------
    # Picard
    RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar -O /opt/picard.jar
    ENV PICARD="/opt/picard.jar"
    
    # GATK
    RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip -O /opt/gatk.zip \
        && unzip /opt/gatk.zip -d /opt/ && rm /opt/gatk.zip
    ENV GATK_HOME=/opt/gatk-4.3.0.0
    ENV PATH="${GATK_HOME}:${PATH}"
    
    # SnpEff
    RUN wget -q http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O /opt/snpeff.zip \
        && unzip /opt/snpeff.zip -d /opt/ && rm /opt/snpeff.zip
    ENV SNPEFF_HOME=/opt/snpEff
    ENV PATH="${SNPEFF_HOME}:${PATH}"
    
    # ANNOVAR
    COPY bin/annovar.latest.tar /opt/
    RUN tar -xvf /opt/annovar.latest.tar -C /opt/ && rm /opt/annovar.latest.tar
    ENV ANNOVAR_HOME=/opt/annovar
    ENV PATH="${ANNOVAR_HOME}:${PATH}"
    
    # -----------------------------
    # Django static files & migrations
    # -----------------------------
    RUN mkdir -p /app/staticfiles && chmod -R 755 /app/staticfiles
    RUN python manage.py collectstatic --noinput
    RUN python manage.py migrate --noinput
    
    # -----------------------------
    # Expose port and run server
    # -----------------------------
    EXPOSE 8000
    CMD ["gunicorn", "omnibioai.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "3"]
    