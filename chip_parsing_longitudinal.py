import math
import csv
import json
import gzip
import argparse
import os
import statistics
import sys
from functools import reduce
import chip_cfg
import csv

from TargetPanel import TargetPanel
from TargetPanelMySQLConnector import TargetPanelMySQLConnector
def parse(config, annotated_vcf, prefix_outfile, translation_filename, bypass_vaf_filter=False):
    def get_longitudinal_pairs(filename):
        sample_2_seqn = dict()
        seqn_2_samples = dict()
        try:
            with open(filename, 'r') as data:
                for file_line in data:
                    seqn, *visitas = file_line.strip().split("\t")
                    sample_2_seqn.update({sampleid: seqn for sampleid in visitas})
                    seqn_2_samples.update({seqn: tuple(visitas)})

            return sample_2_seqn, seqn_2_samples
        except OSError as e:
            print(f"Unable to open {filename}: {e}", file=sys.stderr)
            exit(1)

    def is_in_vaf_threshold_exception_gene():
        return anotacion.get('SYMBOL') in VAF_THRESHOLD_EXCEPTION_GENES

    def get_sample_data(all_alleles, sample_names, vcf_format, vcf_samples):
        def get_sample_info(vcf_format, sdata):
            def extract_gt():
                gt_field = sdata_named.get('GT')
                gt_alleles_idx = [int(gt_field[idx]) for idx in range(len(gt_field)) if idx % 2 == 0]
                gt_alleles = [all_alleles[idx] for idx in gt_alleles_idx]
                return gt_alleles

            def extract_ad():
                res = dict()
                gt_field = sdata_named.get('AD').split(',')
                i = 0
                for allele in all_alleles:  # La AF se anota para los REF + ALT
                    value = int(gt_field[i]) if not gt_field[i] in ['.'] else 0
                    res.update({allele: {'AD': value}})
                    i += 1
                return res

            def extract_af():
                res = dict()
                gt_field = sdata_named.get('AF').split(',')
                af_aggregation = 0.0  # Add up AF of alternate alleles to calculate REF AF
                i = 0
                for allele in all_alleles[1:]:  # La AF se anota solo para los ALT
                    value = float(gt_field[i]) if not gt_field[i] in ['.'] else 0
                    res.update({allele: {'AF': value}})
                    af_aggregation += value
                    i += 1
                # Append REF AF
                res.update({all_alleles[0]: {'AF': 1 - af_aggregation}})

                return res

            def extract_f1r2():
                res = dict()
                gt_field = sdata_named.get('F1R2').split(',')
                i = 0
                for allele in all_alleles:  # Se anota para los REF + ALT
                    value = int(gt_field[i]) if not gt_field[i] in ['.'] else 0
                    res.update({allele: {'F1R2': value}})
                    i += 1
                return res

            def extract_f2r1():
                res = dict()
                gt_field = sdata_named.get('F2R1').split(',')
                i = 0
                for allele in all_alleles:  # Se anota para los REF + ALT
                    value = int(gt_field[i]) if not gt_field[i] in ['.'] else 0
                    res.update({allele: {'F2R1': value}})
                    i += 1
                return res

            def extract_sb():
                res = dict()
                gt_field = sdata_named.get('SB').split(',')
                i = 0
                for allele in all_alleles:  # Se anota para los REF + ALT
                    item = next(iter(gt_field), '.')
                    try:
                        value = int(item) if not item in ['.'] else 0
                    except IndexError:
                        print("aqui")
                    res.update({allele: {'SB': value}})
                    i += 1
                return res

            sdata_named = dict(zip(vcf_format.split(':'), sdata.split(':')))

            res = dict()
            if not sdata_named.get('GT').startswith('.'):
                # Hay llamada, no es NOCALL
                for allele in all_alleles:
                    res.update({allele:
                                    {**extract_af().get(allele, 0),
                                     **extract_ad().get(allele, 0),
                                     **extract_f1r2().get(allele, 0),
                                     **extract_f2r1().get(allele, 0),
                                     **extract_sb().get(allele, 0)}
                                })
            return res

        return {sample: k for i, sample in enumerate(sample_names) if
                (k := get_sample_info(vcf_format, vcf_samples[i]))}

    def extract_CSQ_header(line):
        idx_start = line.find('Allele')
        return line[idx_start:-2].split('|')

    def get_sample_name(samples_vcf, translation_filename=None):
        def create_translation_table_from_file(translation_filename):
            res = dict()
            try:
                with open(translation_filename, 'r') as translation_file:
                    for line in translation_file:
                        orig, dest = line.strip().split('\t')
                        res[orig] = dest
                return res
            except OSError as e:
                print(f"Unable to open {translation_filename}: {e}", file=sys.stderr)
                exit(1)

        if translation_filename:
            translation_table = create_translation_table_from_file(translation_filename)
        else:
            translation_table = dict()

        return [translation_table.get(vcf_id, vcf_id) for vcf_id in samples_vcf]

    def get_whitelist_ncl():
        query = "SELECT CONCAT_WS('_', CHR, POS, REF, VAR) FROM WHITELIST_NUCLEOTIDE"
        return [record[0] for record in chip_cfg.exec_query(db, query)]

    def is_whitelist_ncl(varid):
        return var_id in whitelist_ncl

    def get_whitelist_aa():
        query = "SELECT FEATURE_PREFIX, AA FROM WHITELIST_AMINOACID"
        return {record[0]: int(record[1]) for record in chip_cfg.exec_query(db, query)}

    def is_whitelist_aa():

        def get_aa_position():
            return anotacion.get('Protein_position', '').split('-')[0]

        feature_prefix = anotacion.get('Feature', '').split('.')[0]
        return whitelist_aa.get(feature_prefix) == get_aa_position()

    def get_cosmic():
        res = dict()
        try:
            file = config.get('DB', 'file_COSMIC')
            for line in open(file):
                fields = line.split()
                if len(fields) == 4:
                    table_key = f"{fields[0]}_{fields[1]}_{fields[2]}"
                    res[table_key] = int(fields[3])
            return res
        except Exception as e:
            print(f'ERROR loading COSMIC table : {e}')
            exit(1)

    def in_cosmic(ids, ref, var, cosmic_table):

        def check_cosmic(cosmic_id):
            return cosmic_table.get(cosmic_id).get(ref).get(var) is not None

        return reduce(lambda x, y: x or y, [check_cosmic(cosmic_id) for cosmic_id in ids])

    def get_canonical_transcripts():
        query = "SELECT distinct(FEATURE) from TARGET_PANEL"
        return [record[0] for record in chip_cfg.exec_query(db, query)]

    def is_deleterious(a):

        def has_impact(annotation):
            impact = annotation['IMPACT']

            if not impact:
                print('ERROR: No hay campo IMPACT')
                exit(1)
            else:
                res = impact in ('HIGH', 'MODERATE')

            return res

        return has_impact(a) and a.get('SYMBOL') in panel_genes

    def is_synonymous(a):
        return a['IMPACT'] == 'LOW'

    def get_annotation():
        def get_var_data():
            def unzip_info():
                def unzip_single_info_item(item):
                    ''' Devuelve un dict para un item del campo INFO '''

                    tmp = item.split('=', 1)
                    if len(tmp) == 1:
                        res = {tmp[0]: True}
                    else:
                        # Si es el campo CSQ, lo desglosamos
                        if tmp[0] == 'CSQ':
                            # Obtenemos las anotaciones de los diferentes transcritos por cada alelo
                            csq = dict()
                            for tmp_item_annotation in tmp[1].split(','):
                                item_annotation = dict(zip(csq_header, tmp_item_annotation.split('|')))
                                allele = item_annotation.get('Allele')
                                if item_annotation.get('Feature_type') == 'Transcript':
                                    csq.setdefault(allele, []).append(item_annotation)
                                else:
                                    csq.setdefault(allele, [])
                            res = {tmp[0]: csq}
                        else:
                            res = {tmp[0]: tmp[1]}

                    return res

                res = dict()
                {res.update(unzip_single_info_item(item)) for item in vcf_info.split(';')}
                return res

            def match_csq_allele_annotation(csq):
                def is_padded():
                    return all(alt_allele[idx_padding] == vcf_ref[idx_padding] for alt_allele in vcf_alt.split(','))

                def update_allele(allele):
                    padding_char = vcf_ref[0]
                    if vcf_ref in ['<', '>'] or vcf_alt in ['<', '>']:
                        print("ERROR: Symbolic Allele")
                        exit(1)
                    elif allele == '-':
                        new_alt = padding_char
                    else:
                        new_alt = padding_char + allele

                    if vcf_alt != new_alt:
                        print(f'ERROR: CSQ annotation do not match VCF ALT field : '
                              f'probably VCF normalized after annotation.')
                        exit(1)
                    else:
                        return new_alt

                if (len(vcf_ref) > len(vcf_alt)) and (vcf_ref[0] != vcf_alt[0]):
                    print("Deleccion sin padding")
                    exit(1)
                elif (len(vcf_ref) < len(vcf_alt)) and (vcf_ref[0] != vcf_alt[0]):
                    print("Insercion sin padding")
                    exit(1)
                # elif len(vcf_ref) > 1 and len(vcf_alt) > 1:
                #    print("Complex substution")

                res = dict()
                idx_padding = -1 if vcf_pos == 1 else 0

                if idx_padding == -1:
                    print('ERROR: Padding en posicion 1 del chromosoma')
                    exit(1)
                if is_padded():
                    res = {update_allele(allele): annotation for allele, annotation in csq.items()}
                else:
                    res = csq

                return res

            annotations = unzip_info()
            annotations['CSQ'] = match_csq_allele_annotation(annotations.get('CSQ'))

            return annotations

        def get_canonical():
            def is_canonical(a):
                feature_prefix = a.get('Feature').split('.')[0]
                return feature_prefix in canonical_transcripts

            canonical = None
            for a in anotaciones:
                if is_canonical(a) and a.get('SYMBOL') in panel_genes and is_more_deleterious(a, canonical):
                    canonical = a
            return canonical

        def get_more_deleterious():
            anotacion = None
            for a in anotaciones:
                if a.get('SYMBOL') in panel_genes and is_more_deleterious(a, anotacion):
                    anotacion = a
            return anotacion

        def is_more_deleterious(x_annot, y_annot):
            consequences = {
                'HIGH': 40,
                'MODERATE': 30,
                'LOW': 20,
                'MODIFIER': 10
            }

            x = 0 if x_annot is None else consequences.get(x_annot.get('IMPACT', None), 0)
            y = 0 if y_annot is None else consequences.get(y_annot.get('IMPACT', None), 0)

            return x > y

        nonlocal canonical_transcripts
        var_info = get_var_data()
        anotaciones = var_info['CSQ'].get(variante, 'not_annotated')

        if anotaciones == 'no_annotated':
            print(f'ERROR: No hay ninguna anotacion para el alelo especificado: {var_id}')
            exit(1)
        elif not anotaciones:
            return None
        else:
            canonical = get_canonical()
            if canonical:
                anotacion = canonical
            else:
                anotacion = get_more_deleterious()

            if anotacion:
                anotacion['OTHER_DELETERIOUS'] = [anot.get('HGVSp') for anot in filter(is_deleterious, anotaciones)
                                                  if anot.get('HGVSp', None)]
                anotacion['OTHER_SYNONYMOUS'] = [anot.get('HGVSp') for anot in filter(is_synonymous, anotaciones)
                                                 if anot.get('HGVSp', None)]

            return anotacion

    def get_mutation(alt_field):
        item = alt_field.split(',')

        if len(item) > 1:
            print("ERROR: VCF not normalized. More than 1 alt alleles in the VCF")
            exit(1)
        else:
            return item[0]

    def get_samples_with_mutation(samples):
        return [(sample_id, elem_vaf) for sample_id, sample in samples.items()
                if ((elem_vaf := sample.get(variante).get('AF', 0)) > 0)]

    def get_samples_with_vaf_like_germinal(samples):
        return [id for id, vaf in samples if is_vaf_germ(vaf, nocall_as_germ=False)]

    def heuristic_filter(cosmic_match):

        def do_filter():
            def mean_vaf_cond():
                vaf_samples = [vaf for sample_id, vaf in samples_with_mutation]
                return statistics.mean(vaf_samples) < float(config.get('PARSING', 'HEURISTIC_VAF_THRESHOLD'))

            def num_threshold_cond():
                num_samples_threshold = math.ceil(num_samples
                                                  * float(config.get('PARSING', 'HEURISTIC_FRACTION_SAMPLES')))
                return len(samples_with_mutation) > num_samples_threshold

            return mean_vaf_cond() and num_threshold_cond()

        def common_in_cosmic():
            def test(num_matches):
                num_matches > int(config.get('PARSING', 'HEURISTIC_COSMIC_OCURRENCES_THRESHOLD'))

            return reduce(lambda x, y: x or y, [test(i[1]) for i in cosmic_match]) if cosmic_match else False

        return do_filter() and not common_in_cosmic()

    def get_previously_identified():
        query = "SELECT CONCAT_WS('_', CHR, POS, REF, VAR) FROM PREVIOUSLY_IDENTIFIED"
        return [record[0] for record in chip_cfg.exec_query(db, query)]

    def get_artifacts():
        query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM INTERNALLY_IDENTIFIED_NO_DRIVER GROUP BY VAR_ID"
        return {record[0]: record[1] for record in chip_cfg.exec_query(db, query)}

    def is_artifact():
        return var_id in artifacts

    def is_previously_identified():
        return var_id in previously_identified

    def get_internally_identified():
        query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM INTERNALLY_IDENTIFIED GROUP BY VAR_ID"
        return {record[0]: record[1] for record in chip_cfg.exec_query(db, query)}

    def create_printing_deleterious():
        def get_cosmic_match(ids):
            return [(cosmic_id, cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}"))
                    for cosmic_id in ids if cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}")]

        def predictor_splitting(prediction_chunk):
            tmp = prediction_chunk.split('(')

            if prediction_chunk:
                try:
                    desc = tmp[0]
                    score = tmp[1].split(')')[0]
                except IndexError:
                    desc = prediction_chunk
                    score = ''
            else:
                desc = ''
                score = ''
            return desc, score

        def get_existing_variation(anotacion):
            def get_ids(id_chunk, mask_chunk):
                # mask = iter(mask_chunk.split('&'))
                mask = mask_chunk.split('&')
                ids = id_chunk.split('&')

                return [ident for ident, flag in zip(ids, mask) if flag == '1']

            somatic = get_ids(anotacion.get('Existing_variation'), anotacion.get('SOMATIC'))
            pheno = get_ids(anotacion.get('Existing_variation'), anotacion.get('PHENO'))

            return somatic, pheno

        sift_desc, sift_score = predictor_splitting(anotacion.get('SIFT'))
        polyphen_desc, polyphen_score = predictor_splitting(anotacion.get('PolyPhen'))
        somatic_ids, pheno_ids = get_existing_variation(anotacion)
        cosmic_match = get_cosmic_match(somatic_ids)

        item_var = {
            'FILTER': vcf_filter,
            'VAR_ID': var_id,
            'CHR': vcf_chrom,
            'POS': vcf_pos,
            'REF': vcf_ref,
            'ALT': variante
        }

        item_samples = {
            'LONGITUDINAL_SEQN': seqn,
            'SAMPLE_VIS1': vis1,
            'REF_DEPTH_VIS1': samples_gt.get(vis1, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('AD'),
            'F1R2_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F2R1'),
            'SAMPLE_VIS3': vis3,
            'REF_DEPTH_VIS3': samples_gt.get(vis3, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('AD'),
            'F1R2_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F2R1'),
        }

        item_vaf = {
            'VAF_VIS1': vaf_vis1 if vaf_vis1 is not None else 'NOCALL',
            'VAF_VIS3': vaf_vis3 if vaf_vis3 is not None else 'NOCALL',
            'VAF_ABSOLUTA': (vaf_vis3 - vaf_vis1) if vaf_vis1 is not None and vaf_vis3 is not None else 'NA',
            'VAF_RELATIVA': (((vaf_vis3 - vaf_vis1) / vaf_vis1) * 100)
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA',
            'VAF_RATIO': vaf_vis3 / vaf_vis1
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA'
        }

        item_anot = {
            'GENE': anotacion.get('SYMBOL'),
            'ALLELE': anotacion.get('Allele'),
            'VARIANT_CLASS': anotacion.get('VARIANT_CLASS'),
            'BIOTYPE': anotacion.get('BIOTYPE'),
            'TYPE': anotacion.get('Feature_type'),
            'FEATURE': anotacion.get('Feature'),
            'EXON': anotacion.get('EXON'),
            'HGVSc': anotacion.get('HGVSc'),
            'HGVSp': anotacion.get('HGVSp'),
            'cDNA_position': anotacion.get('cDNA_position'),
            'CDS_position': anotacion.get('CDS_position'),
            'Protein_position': anotacion.get('Protein_position'),
            'Amino_acids': anotacion.get('Amino_acids'),
            'Codons': anotacion.get('Codons'),
            'STRAND': anotacion.get('STRAND'),
            'REFSEQ_MATCH': anotacion.get('REFSEQ_MATCH'),
            'GIVEN_REF': anotacion.get('GIVEN_REF'),
            'USED_REF': anotacion.get('USED_REF'),
            'BAM_EDIT': anotacion.get('BAM_EDIT'),
            'IMPACT': anotacion.get('IMPACT'),
            'CONSEQUENCE': anotacion.get('Consequence'),
            'SOMATIC': ';'.join(somatic_ids),
            'PHENO': ';'.join(pheno_ids),
            'CLIN_SIG': anotacion.get('CLIN_SIG'),
            'COSMIC_MATCH': ';'.join([i[0] for i in cosmic_match]),
            'HEURISTIC': 'Y' if heuristic_filter(cosmic_match) else 'N',
            'PREVIOUSLY_IDENTIFIED': 'Y' if is_previously_identified() else 'N',
            'WHITELIST': 'Y' if is_whitelisted else 'N',
            'INTERNALLY_IDENTIFIED': internally_identified.get(var_id),
            'ARTIFACT': 'Y' if is_artifact() else 'N'
        }

        item_predictors = {
            'SIFT_DESC': sift_desc,
            'SIFT_SCORE': sift_score,
            'POLYPHEN_DESC': polyphen_desc,
            'POLYPHEN_SCORE': polyphen_score,
            'CADD_PHRED': anotacion.get('CADD_phred'),
            'CADD_RAW': anotacion.get('CADD_raw'),
        }
        return {**item_var, **item_anot, **item_predictors, **item_samples, **item_vaf}

    def create_printing_synonymous():

        def get_cosmic_match(ids):
            return [(cosmic_id, cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}"))
                    for cosmic_id in ids if cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}")]

        def get_existing_variation(anotacion):
            """
                Obtiene los IDs en las BD externas
            """

            def get_ids(id_chunk, mask_chunk):
                # mask = iter(mask_chunk.split('&'))
                mask = mask_chunk.split('&')
                ids = id_chunk.split('&')

                return [ident for ident, flag in zip(ids, mask) if flag == '1']

            somatic = get_ids(anotacion.get('Existing_variation'), anotacion.get('SOMATIC'))
            pheno = get_ids(anotacion.get('Existing_variation'), anotacion.get('PHENO'))

            return somatic, pheno

        somatic_ids, pheno_ids = get_existing_variation(anotacion)
        cosmic_match = get_cosmic_match(somatic_ids)

        item_var = {
            'FILTER': vcf_filter,
            'VAR_ID': var_id,
            'CHR': vcf_chrom,
            'POS': vcf_pos,
            'REF': vcf_ref,
            'ALT': variante
        }

        item_samples = {
            'LONGITUDINAL_SEQN': seqn,
            'SAMPLE_VIS1': vis1,
            'REF_DEPTH_VIS1': samples_gt.get(vis1, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('AD'),
            'F1R2_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F2R1'),
            'SAMPLE_VIS3': vis3,
            'REF_DEPTH_VIS3': samples_gt.get(vis3, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('AD'),
            'F1R2_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F2R1'),
        }

        item_vaf = {
            'VAF_VIS1': vaf_vis1 if vaf_vis1 is not None else 'NOCALL',
            'VAF_VIS3': vaf_vis3 if vaf_vis3 is not None else 'NOCALL',
            'VAF_ABSOLUTA': (vaf_vis3 - vaf_vis1) if vaf_vis1 is not None and vaf_vis3 is not None else 'NA',
            'VAF_RELATIVA': (((vaf_vis3 - vaf_vis1) / vaf_vis1) * 100)
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA',
            'VAF_RATIO': vaf_vis3 / vaf_vis1
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA'
        }

        item_anot = {
            'GENE': anotacion.get('SYMBOL'),
            'ALLELE': anotacion.get('Allele'),
            'VARIANT_CLASS': anotacion.get('VARIANT_CLASS'),
            'BIOTYPE': anotacion.get('BIOTYPE'),
            'TYPE': anotacion.get('Feature_type'),
            'FEATURE': anotacion.get('Feature'),
            'EXON': anotacion.get('EXON'),
            'HGVSc': anotacion.get('HGVSc'),
            'HGVSp': anotacion.get('HGVSp'),
            'cDNA_position': anotacion.get('cDNA_position'),
            'CDS_position': anotacion.get('CDS_position'),
            'Protein_position': anotacion.get('Protein_position'),
            'Amino_acids': anotacion.get('Amino_acids'),
            'Codons': anotacion.get('Codons'),
            'STRAND': anotacion.get('STRAND'),
            'REFSEQ_MATCH': anotacion.get('REFSEQ_MATCH'),
            'GIVEN_REF': anotacion.get('GIVEN_REF'),
            'USED_REF': anotacion.get('USED_REF'),
            'BAM_EDIT': anotacion.get('BAM_EDIT'),
            'IMPACT': anotacion.get('IMPACT'),
            'CONSEQUENCE': anotacion.get('Consequence'),
            'SOMATIC': ';'.join(somatic_ids),
            'PHENO': ';'.join(pheno_ids),
            'CLIN_SIG': anotacion.get('CLIN_SIG'),
            'COSMIC_MATCH': ';'.join([i[0] for i in cosmic_match]),
            'HEURISTIC': 'Y' if heuristic_filter(cosmic_match) else 'N',
            'PREVIOUSLY_IDENTIFIED': 'Y' if is_previously_identified() else 'N',
            'WHITELIST': 'Y' if is_whitelisted else 'N',
            'INTERNALLY_IDENTIFIED': internally_identified.get(var_id),
            'ARTIFACT': 'Y' if is_artifact() else 'N'
        }

        return {**item_var, **item_anot, **item_samples, **item_vaf}

    def create_printing_dump():

        def get_cosmic_match(ids):
            return [(cosmic_id, cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}"))
                    for cosmic_id in ids if cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}")]

        def predictor_splitting(prediction_chunk):
            tmp = prediction_chunk.split('(')

            if prediction_chunk:
                try:
                    desc = tmp[0]
                    score = tmp[1].split(')')[0]
                except IndexError:
                    desc = prediction_chunk
                    score = ''
            else:
                desc = ''
                score = ''
            return desc, score

        def get_existing_variation(anotacion):
            def get_ids(id_chunk, mask_chunk):
                # mask = iter(mask_chunk.split('&'))
                mask = mask_chunk.split('&')
                ids = id_chunk.split('&')

                return [ident for ident, flag in zip(ids, mask) if flag == '1']

            somatic = get_ids(anotacion.get('Existing_variation'), anotacion.get('SOMATIC'))
            pheno = get_ids(anotacion.get('Existing_variation'), anotacion.get('PHENO'))

            return somatic, pheno

        sift_desc, sift_score = predictor_splitting(anotacion.get('SIFT'))
        polyphen_desc, polyphen_score = predictor_splitting(anotacion.get('PolyPhen'))
        somatic_ids, pheno_ids = get_existing_variation(anotacion)
        cosmic_match = get_cosmic_match(somatic_ids)

        item_var = {
            'FILTER': vcf_filter,
            'VAR_ID': var_id,
            'CHR': vcf_chrom,
            'POS': vcf_pos,
            'REF': vcf_ref,
            'ALT': variante
        }

        item_samples = {
            'LONGITUDINAL_SEQN': seqn,
            'SAMPLE_VIS1': vis1,
            'REF_DEPTH_VIS1': samples_gt.get(vis1, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('AD'),
            'F1R2_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS1': samples_gt.get(vis1, {}).get(variante, {}).get('F2R1'),
            'SAMPLE_VIS3': vis3,
            'REF_DEPTH_VIS3': samples_gt.get(vis3, {}).get(vcf_ref, {}).get('AD'),
            'ALT_DEPTH_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('AD'),
            'F1R2_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F1R2'),
            'F2R1_VIS3': samples_gt.get(vis3, {}).get(variante, {}).get('F2R1'),
        }

        item_vaf = {
            'VAF_VIS1': vaf_vis1 if vaf_vis1 is not None else 'NOCALL',
            'VAF_VIS3': vaf_vis3 if vaf_vis3 is not None else 'NOCALL',
            'VAF_ABSOLUTA': (vaf_vis3 - vaf_vis1) if vaf_vis1 is not None and vaf_vis3 is not None else 'NA',
            'VAF_RELATIVA': (((vaf_vis3 - vaf_vis1) / vaf_vis1) * 100)
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA',
            'VAF_RATIO': vaf_vis3 / vaf_vis1
            if vaf_vis1 is not None and vaf_vis3 is not None and vaf_vis1 > 0 else 'NA'
        }

        item_anot = {
            'GENE': anotacion.get('SYMBOL'),
            'ALLELE': anotacion.get('Allele'),
            'VARIANT_CLASS': anotacion.get('VARIANT_CLASS'),
            'BIOTYPE': anotacion.get('BIOTYPE'),
            'TYPE': anotacion.get('Feature_type'),
            'FEATURE': anotacion.get('Feature'),
            'EXON': anotacion.get('EXON'),
            'HGVSc': anotacion.get('HGVSc'),
            'HGVSp': anotacion.get('HGVSp'),
            'cDNA_position': anotacion.get('cDNA_position'),
            'CDS_position': anotacion.get('CDS_position'),
            'Protein_position': anotacion.get('Protein_position'),
            'Amino_acids': anotacion.get('Amino_acids'),
            'Codons': anotacion.get('Codons'),
            'STRAND': anotacion.get('STRAND'),
            'REFSEQ_MATCH': anotacion.get('REFSEQ_MATCH'),
            'GIVEN_REF': anotacion.get('GIVEN_REF'),
            'USED_REF': anotacion.get('USED_REF'),
            'BAM_EDIT': anotacion.get('BAM_EDIT'),
            'IMPACT': anotacion.get('IMPACT'),
            'CONSEQUENCE': anotacion.get('Consequence'),
            'SOMATIC': ';'.join(somatic_ids),
            'PHENO': ';'.join(pheno_ids),
            'CLIN_SIG': anotacion.get('CLIN_SIG'),
            'COSMIC_MATCH': ';'.join([i[0] for i in cosmic_match]),
            'HEURISTIC': 'Y' if heuristic_filter(cosmic_match) else 'N',
            'PREVIOUSLY_IDENTIFIED': 'Y' if is_previously_identified() else 'N',
            'WHITELIST': 'Y' if is_whitelisted else 'N',
            'INTERNALLY_IDENTIFIED': internally_identified.get(var_id),
            'ARTIFACT': 'Y' if is_artifact() else 'N'
        }

        item_predictors = {
            'SIFT_DESC': sift_desc,
            'SIFT_SCORE': sift_score,
            'POLYPHEN_DESC': polyphen_desc,
            'POLYPHEN_SCORE': polyphen_score,
            'CADD_PHRED': anotacion.get('CADD_phred'),
            'CADD_RAW': anotacion.get('CADD_raw'),
        }
        return {**item_var, **item_anot, **item_predictors, **item_samples, **item_vaf}

    def is_vaf_germ(vaf, nocall_as_germ=False):
        if vaf is None or vaf == '.':
            res = True if nocall_as_germ else False
        else:
            vaf_casted = float(vaf)
            res = (0.45 < vaf_casted < 0.55) or vaf_casted > 0.85

        return res

    def seems_somatic():
        def is_in_population():
            max_af = anotacion.get('MAX_AF', -1)
            if max_af == -1:
                print('ERROR: VEP was executed without --max_af option')
                exit(1)
            elif not max_af:
                return False
            else:
                return float(max_af) > MAX_MAF_THRESHOLD

        ok_samples_with_mutation = len(samples_with_mutation) < MAX_NUM_SAMPLES_WITH_VAR
        ok_samples_with_vaf_like_germinal = len(samples_with_vaf_like_germinal) < MAX_NUM_SAMPLES_GERMINAL_VAF

        return ok_samples_with_mutation and ok_samples_with_vaf_like_germinal and not is_in_population()

    longitudinal_distribution_dump = config.getboolean('GENERAL', 'longitudinal_distribution_dump')

    bypass_cosmic = False

    VAF_THRESHOLD = config.getfloat('PARSING', 'VAF_THRESHOLD')
    MAX_MAF_THRESHOLD = config.getfloat('PARSING', 'MAX_MAF_THRESHOLD')
    MAX_NUM_SAMPLES_GERMINAL_VAF = config.getint('PARSING', 'MAX_NUM_SAMPLES_GERMINAL_VAF')
    VAF_THRESHOLD_EXCEPTION_GENES = config.get('PARSING', 'VAF_THRESHOLD_EXCEPTION_GENES').split(',')

    db = chip_cfg.connect_db(config)

    path_outfile = os.path.split(prefix_outfile)[0]
    os.makedirs(path_outfile, exist_ok=True)
    # Output files
    header_flag_deleterious = True
    out_deletereous = open(f'{prefix_outfile}.candidates.tsv', 'wt')
    tsv_deletereous = csv.writer(out_deletereous, delimiter='\t')

    header_flag_synonymous = True
    out_synonymous = open(f'{prefix_outfile}.synonymous.tsv', 'wt')
    tsv_synonymous = csv.writer(out_synonymous, delimiter='\t')

    if longitudinal_distribution_dump:
        header_flag_dump = True
        out_dump = open(f'{prefix_outfile}.dump.tsv', 'wt')
        tsv_dump = csv.writer(out_dump, delimiter='\t')

    longitudinal_pairs = config.get('DIR', 'longitudinal_pairs')
    sample_2_seqn, seqn_2_samples = get_longitudinal_pairs(longitudinal_pairs)

    cosmic_table = get_cosmic() if not bypass_cosmic else {}
    whitelist_ncl = get_whitelist_ncl()
    whitelist_aa = get_whitelist_aa()
    previously_identified = get_previously_identified()
    internally_identified = get_internally_identified()
    artifacts = get_artifacts()

    # Lee datos del panel
    db_connection = TargetPanelMySQLConnector(config)
    panel = TargetPanel(config, db_connection)
    panel_genes = panel.genes()
    canonical_transcripts = panel.transcripts()

    with gzip.open(annotated_vcf, 'rt') as vcf_file:
        for line_chunk in vcf_file:
            if line_chunk.startswith('##INFO=<ID=CSQ'):
                # Extract CSQ header from vcf
                csq_header = extract_CSQ_header(line_chunk)
            elif line_chunk.startswith('#CHROM'):
                sample_names = get_sample_name(line_chunk.split()[9:], translation_filename)
                num_samples = len(sample_names)
                MAX_NUM_SAMPLES_WITH_VAR = math.ceil(
                    num_samples * float(config.get('PARSING', 'TOTAL_SAMPLES_FRACTION')))
            elif not line_chunk.startswith('#'):
                vcf_chrom, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, *vcf_samples \
                    = line_chunk.strip().split('\t')
    
                variante = get_mutation(vcf_alt)
                var_id = f"""{vcf_chrom}_{vcf_pos}_{vcf_ref}_{variante}"""

                anotacion = get_annotation()

                if anotacion:
                    samples_gt = get_sample_data([vcf_ref, *vcf_alt.split(',')], sample_names, vcf_format,
                                                 vcf_samples)

                    samples_with_mutation = get_samples_with_mutation(samples_gt)

                    samples_with_vaf_like_germinal = get_samples_with_vaf_like_germinal(samples_with_mutation)

                    is_whitelisted = is_whitelist_ncl(var_id) or is_whitelist_aa()
                    if is_whitelisted or seems_somatic():
                        for seqn in set([seqn for sampleid, vaf in samples_with_mutation
                                         if (seqn := sample_2_seqn.get(sampleid, False))]):
                            vis1, vis3 = seqn_2_samples.get(seqn)

                            if ((tmp_vis1 := samples_gt.get(vis1, {}).get(variante, False))
                                    and (tmp_vis3 := samples_gt.get(vis3, {}).get(variante, False))):
                                vaf_vis1 = float(tmp_vis1.get('AF')) if tmp_vis1 else None
                                vaf_vis3 = float(tmp_vis3.get('AF')) if tmp_vis3 else None

                                if not is_vaf_germ(vaf_vis1, nocall_as_germ=True) or \
                                        not is_vaf_germ(vaf_vis3, nocall_as_germ=True):
                                    if longitudinal_distribution_dump:
                                        item = create_printing_dump()
                                        if header_flag_dump:
                                            tsv_dump.writerow(item.keys())
                                            header_flag_dump = False
                                        tsv_dump.writerow(item.values())

                                    if not bypass_vaf_filter:
                                        ok_vaf1 = vaf_vis1 >= VAF_THRESHOLD if vaf_vis1 is not None else False
                                        ok_vaf3 = vaf_vis3 >= VAF_THRESHOLD if vaf_vis3 is not None else False

                                        ok_vaf = ok_vaf1 or ok_vaf3 or is_in_vaf_threshold_exception_gene()
                                    else:
                                        ok_vaf = True

                                    if ok_vaf:
                                        if is_deleterious(anotacion):
                                            item = create_printing_deleterious()
                                            if header_flag_deleterious:
                                                tsv_deletereous.writerow(item.keys())
                                                header_flag_deleterious = False
                                            tsv_deletereous.writerow(item.values())
                                        elif is_synonymous(anotacion):
                                            item = create_printing_synonymous()
                                            if header_flag_synonymous:
                                                tsv_synonymous.writerow(item.keys())
                                                header_flag_synonymous = False
                                            tsv_synonymous.writerow(item.values())
    out_deletereous.close()
    out_synonymous.close()
    db.close()
