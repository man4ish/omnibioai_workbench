import os
import subprocess
from collections import Counter

from django.shortcuts import render, redirect
from django.contrib import messages
from django.conf import settings

# Set paths (adjust as per your environment)
UPLOAD_FOLDER = os.path.join(settings.BASE_DIR, 'data')
ANNOVAR_PATH = '/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/annovar'
HUMANDB_PATH = os.path.join(ANNOVAR_PATH, 'humandb')

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)


def summarize_variants_by_chromosome(variants):
    chroms = [v.get('Chr', 'Unknown') for v in variants]
    counts = Counter(chroms)
    return dict(counts)


class BaseAnnotator:
    def __init__(self, vcf_file_path, output_prefix):
        self.vcf_file_path = vcf_file_path
        self.output_prefix = output_prefix

    def run(self):
        raise NotImplementedError

    def parse_output(self):
        raise NotImplementedError


class AnnovarAnnotator(BaseAnnotator):
    def run(self):
        cmd = [
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'),
            self.vcf_file_path,
            HUMANDB_PATH,
            '-buildver', 'hg19',
            '-out', self.output_prefix,
            '-remove',
            '-protocol', 'refGene,clinvar_20210501',
            '-operation', 'g,f',
            '-nastring', '.',
            '-vcfinput'
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"ANNOVAR annotation failed: {result.stderr}")

    def parse_output(self):
        variants = []
        filepath = f"{self.output_prefix}.hg19_multianno.txt"
        if not os.path.exists(filepath):
            return variants

        with open(filepath, 'r') as f:
            headers = f.readline().strip().split('\t')
            for line in f:
                fields = line.strip().split('\t')
                variant = dict(zip(headers, fields))
                variants.append(variant)
        return variants


class VepAnnotator(BaseAnnotator):
    def run(self):
        self.annotated_vcf = f"{self.output_prefix}.vep_annotated.vcf"

        vep_script = "/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/ensembl-vep-release-114/vep"

        if not os.path.exists(vep_script):
            raise FileNotFoundError(f"VEP script not found at: {vep_script}")

        cmd = [
            vep_script,
            "--input_file", self.vcf_file_path,
            "--output_file", self.annotated_vcf,
            "--vcf",
            "--force_overwrite",
            "--offline",
            "--cache",
            "--everything",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise Exception(f"VEP annotation failed:\n{result.stderr}")

    def parse_output(self):
        variants = []
        if not os.path.exists(self.annotated_vcf):
            return variants

        with open(self.annotated_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                info_dict = dict(
                    item.split('=') if '=' in item else (item, '')
                    for item in fields[7].split(';')
                )

                variant = {
                    'Chr': fields[0],
                    'Pos': fields[1],
                    'ID': fields[2],
                    'Ref': fields[3],
                    'Alt': fields[4],
                    'Qual': fields[5],
                    'Filter': fields[6],
                    'Info': info_dict
                }
                variants.append(variant)
        return variants


class SnpEffAnnotator(BaseAnnotator):
    def run(self):
        self.annotated_vcf = f"{self.output_prefix}.snpeff_annotated.vcf"

        if not hasattr(self, 'genome_version'):
            raise AttributeError("Missing 'genome_version'. Set self.genome_version to something like 'GRCh38.99'.")

        snpeff_jar = "/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/snpEff/snpEff.jar"

        if not os.path.exists(snpeff_jar):
            raise FileNotFoundError(f"snpEff jar not found at {snpeff_jar}")

        cmd = [
            "java", "-Xmx4g", "-jar", snpeff_jar,
            self.genome_version,
            self.vcf_file_path
        ]

        with open(self.annotated_vcf, 'w') as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            raise Exception(f"snpEff annotation failed:\n{result.stderr}")

    def parse_output(self):
        variants = []
        if not os.path.exists(self.annotated_vcf):
            return variants

        with open(self.annotated_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                info_dict = dict(
                    item.split('=') if '=' in item else (item, '')
                    for item in fields[7].split(';')
                )

                variant = {
                    'Chr': fields[0],
                    'Pos': fields[1],
                    'ID': fields[2],
                    'Ref': fields[3],
                    'Alt': fields[4],
                    'Qual': fields[5],
                    'Filter': fields[6],
                    'Info': info_dict
                }
                variants.append(variant)
        return variants


class BcftoolsAnnotator(BaseAnnotator):
    def run(self):
        self.annotated_vcf = f"{self.output_prefix}.bcftools_annotated.vcf"

        cmd = [
            "bcftools", "annotate",
            "-o", self.annotated_vcf,
            "-O", "v",
            self.vcf_file_path
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"BCFtools annotation failed: {result.stderr}")

    def parse_output(self):
        variants = []
        if not os.path.exists(self.annotated_vcf):
            return variants

        with open(self.annotated_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                info_dict = dict(item.split('=') if '=' in item else (item, '') for item in fields[7].split(';'))

                variant = {
                    'Chr': fields[0],
                    'Pos': fields[1],
                    'ID': fields[2],
                    'Ref': fields[3],
                    'Alt': fields[4],
                    'Qual': fields[5],
                    'Filter': fields[6],
                    'Info': info_dict
                }
                variants.append(variant)
        return variants


def index(request):
    if request.method == 'POST':
        uploaded_file = request.FILES.get('vcf_file')
        if not uploaded_file:
            messages.error(request, 'No file uploaded')
            return redirect('index')

        annotator_choice = request.POST.get('annotator', 'annovar').lower()

        # Save uploaded file
        filepath = os.path.join(UPLOAD_FOLDER, uploaded_file.name)
        with open(filepath, 'wb+') as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        output_prefix = os.path.join(UPLOAD_FOLDER, 'annotated_output')

        annotator_map = {
            'annovar': AnnovarAnnotator,
            'vep': VepAnnotator,
            'snpeff': SnpEffAnnotator,
            'bcftools': BcftoolsAnnotator
        }

        annotator_class = annotator_map.get(annotator_choice)
        if not annotator_class:
            messages.error(request, f"Unknown annotator selected: {annotator_choice}")
            return redirect('index')

        annotator = annotator_class(filepath, output_prefix)

        if annotator_choice == 'snpeff':
            annotator.genome_version = 'GRCh38.99'  # adjust genome version as needed

        try:
            annotator.run()
            variants = annotator.parse_output()
        except Exception as e:
            messages.error(request, f"Annotation failed: {str(e)}")
            return redirect('index')

        variant_counts = summarize_variants_by_chromosome(variants)

        return render(request, 'variant_annotation/result.html', {
            'variants': variants,
            'variant_counts': variant_counts,
            'annotator': annotator_choice
        })

    return render(request, 'variant_annotation/index.html')
