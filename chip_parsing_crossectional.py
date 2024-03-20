import math
import os
import sqlite3
import csv
import json
import gzip
import argparse
import statistics
from configparser import NoOptionError
from functools import reduce
import chip_cfg
import csv
import sys
from TargetPanel import TargetPanel
from TargetPanelMySQLConnector import TargetPanelMySQLConnector

"""
    Read VCF file and make a list of mutations candidates to be curated manually
"""

def parse(config, annotated_vcf, tsv_outfile, translation_filename, bypass_vaf_filter=False):
    def get_whitelist_ncl():
        query = "SELECT CONCAT_WS('_', CHR, POS, REF, VAR) FROM WHITELIST_NUCLEOTIDE"
        return [record[0] for record in chip_cfg.exec_query(db, query)]

    def get_whitelist_aa():
        query = "SELECT FEATURE_PREFIX, AA FROM WHITELIST_AMINOACID"
        return {record[0]: int(record[1]) for record in chip_cfg.exec_query(db, query)}

    def is_in_vaf_threshold_exception_gene():
        """
            Gene list wich do not apply VAF threshold
        """
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

        res = dict()
        i = 0
        for sample in sample_names:
            k = get_sample_info(vcf_format, vcf_samples[i])
            i += 1
            if k:
                res.update({sample: k})
        return res

    def has_impact(annotation):
        impact = annotation['IMPACT']

        if not impact:
            print('ERROR: No hay campo IMPACT')
            exit(1)
        else:
            res = impact in ('HIGH', 'MODERATE')

        return res

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

    # -----
    # WHITE LIST

    def is_whitelist_ncl(varid):
        return var_id in whitelist_ncl

    def is_whitelist_aa():

        def get_aa_position():
            return anotacion.get('Protein_position', '').split('-')[0]

        feature_prefix = anotacion.get('Feature', '').split('.')[0]
        return whitelist_aa.get(feature_prefix) == get_aa_position()

    def get_cosmic():
        """
            Create COSMIC table from the file
        """
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
        """ Load canonical transcripts from database"""

        query = "SELECT distinct(FEATURE) from TARGET_PANEL"
        return [record[0] for record in chip_cfg.exec_query(db, query)]

    def is_deleterious(a):

        def has_impact(annotation):
            '''
                Comprueba si el impacto es HIGH or MODERATE
            '''

            impact = annotation['IMPACT']

            if not impact:
                print('ERROR: No hay campo IMPACT')
                exit(1)
            else:
                res = impact in ('HIGH', 'MODERATE')

            return res

        return has_impact(a) and a.get('SYMBOL') in panel_genes

    def get_annotation():
        """
            Devuelve la 'anotación principal'. Se considera 'anotacion principal' la correspondiente al transcrito
            canonico y ,en su defecto, a la más deleterea de las anotaciones (una de ellas si hay varias)
            Se incluyen campos adicionales con los id de todas las anotaciones deletereas y sinonimas
        """

        def get_var_data():
            """  Devuelve un dict con los datos del campo INFO. El campo CSQ se modifica para que
                la `key` se corresponda con los alelos del campo ALT
            """

            def unzip_info():
                """
                    Devuelve un dict para todos los items del campo INFO
                """

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
                """
                    Devuelve dict donde la `key` se corresponde con el el campo `alt` del VCF.
                    Modifica las `key` si es necesario, y las deja igual si no lo es

                    NOTA de la especificacion VCF:

                        For simple insertions and deletions in which either the REF or one of the
                        ALT alleles would otherwise be null/empty, the REF and ALT Strings must
                        include the base before the event, unless the event
                        occurs at position 1 on the contig in which case it must include the base after
                        the event; this padding base is not required (although it is permitted)
                        for e.g. complex substitutions or other events where all alleles have at least
                        one base represented in their Strings.
                """

                def is_padded():
                    """
                        Devuelve un booleano indicando si TODOS los alt (el el ref) son padded

                        NOTA: Parece (pero no estoy seguro) que VEP solo anota el alelo trimado (eliminado el padding)
                            si todos los alt tienen padding.
                    """
                    return all(alt_allele[idx_padding] == vcf_ref[idx_padding] for alt_allele in vcf_alt.split(','))

                def update_allele(allele):
                    """
                        Modifica `allele` para que coincida con el campo ALT

                        NOTE: No usamos vcf_alt para generar el alelo a devovel para asegurarnos que
                            lo hacemos correctamnte, de modo que si no hay un match entre el calculado
                            (el que devolvemos) y vcf_alt, salte un error
                    """

                    # En este caso, como todavia no se contempla la causistica del padding en la primero posicion
                    # del cromosoma, el caracter de padding es la primera base de vcf_ref y vcf_alt
                    padding_char = vcf_ref[0]
                    if vcf_ref in ['<', '>'] or vcf_alt in ['<', '>']:
                        # Symbolic allele
                        # Desconozco como aparecería en el CSQ; me imagino que igua
                        # Como este caso no lo controlo. Devuelvo ERROR
                        print("ERROR: Symbolic Allele")
                        exit(1)
                    elif allele == '-':
                        # Si el alelo anotado en CSQ es '-' significa que es una deleccion, y por tanto, se debe
                        # establecer la clave de CSQ con el caracter de padding.
                        new_alt = padding_char
                    else:
                        # Hay que modificar la `key` de las anotaciones CSQ porque tb hay padding
                        new_alt = padding_char + allele

                    if vcf_alt != new_alt:
                        # Esta comprobacion tiene sentido porque el VCF esta normalizado y, por tanto solo hay un
                        # elemento en vcf_alt y tiene que coincidir con el anotado
                        print(f'ERROR: CSQ annotation do not match VCF ALT field : '
                              f'probably VCF normalized after annotation.')
                        exit(1)
                    else:
                        return new_alt

                # Hago unas comprobaciones de seguridad
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
                # Mientras no lo manejo bien (hay que repasar update_allele)
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
            """
                Devuelve la anotacion canonica de entre los genes a considerar. Si hay varias devuelve la deleterea. Si hay varias deletereas devuelve
                la mas deleterea o, en su defecto, una cualquiera.
            """

            def is_canonical(a):
                feature_prefix = a.get('Feature').split('.')[0]
                return feature_prefix in canonical_transcripts

            canonical = None
            for a in anotaciones:
                if is_canonical(a) and a.get('SYMBOL') in panel_genes and is_more_deleterious(a, canonical):
                    canonical = a
            return canonical

        def get_more_deleterious():
            """
                Devuelve la anotacion deleterea con mas impacto de entre los genes a considerr.
                Si hay varias de ellas devuelve una.
            """

            anotacion = None
            for a in anotaciones:
                if a.get('SYMBOL') in panel_genes and is_more_deleterious(a, anotacion):
                    anotacion = a
            return anotacion

        def is_more_deleterious(x_annot, y_annot):
            """
                    return TRUE if x is more deletereous than y. Any value is more deleterious than None
                """
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
            # Las anotaciones no se referian a transcritos (p.e a region intergenica)
            return None
        else:
            # Obtiene la anotacion del transcrito canonico
            canonical = get_canonical()
            if canonical:
                anotacion = canonical
            else:
                # Obtiene la anotacion mas deleterea
                anotacion = get_more_deleterious()

            # Si la mutacion cae en un gen excluido (p.e 'LUC7L2') no hay anotacion
            if anotacion:
                anotacion['OTHER_DELETERIOUS'] = [anot.get('HGVSp') for anot in filter(is_deleterious, anotaciones) if
                                                  anot.get('HGVSp', False)]

            return anotacion

    def get_artifacts():
        query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM INTERNALLY_IDENTIFIED_NO_DRIVER GROUP BY VAR_ID"
        return {record[0]: record[1] for record in chip_cfg.exec_query(db, query)}

    def is_artifact():
        return var_id in artifacts

    def get_mutation(alt_field):
        item = alt_field.split(',')

        if len(item) > 1:
            print("ERROR: VCF not normalized. More than 1 alt alleles in the VCF")
            exit(1)
        else:
            return item[0]

    def get_samples_with_mutation(samples):
        return [(sample_id, sample.get(variante).get('AF')) for sample_id, sample in samples.items()
                if sample.get(variante).get('AF') > 0]

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

    def is_previously_identified():
        return var_id in previously_identified

    def get_internally_identified():
        query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM INTERNALLY_IDENTIFIED GROUP BY VAR_ID"
        return {record[0]: record[1] for record in chip_cfg.exec_query(db, query)}

    def create_printing_item():

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
        sample = samples_gt[sample_id]
        # longitudinal_seqn, longitudinal_visita = longitudinal.get(sample_id, (None, None))
        res = {
            'FILTER': vcf_filter,
            "VAR_ID": var_id,
            'SAMPLE': sample_id,
            'GENE': anotacion.get('SYMBOL'),
            'CHR': vcf_chrom,
            'POS': vcf_pos,
            'REF': vcf_ref,
            'ALT': variante,
            'REF_DEPTH': sample.get(vcf_ref).get('AD'),
            'ALT_DEPTH': sample.get(variante).get('AD'),
            'F1R2': sample.get(variante).get('F1R2'),
            'F2R1': sample.get(variante).get('F2R1'),
            'VAF': vaf,
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
            'SIFT_DESC': sift_desc,
            'SIFT_SCORE': sift_score,
            'POLYPHEN_DESC': polyphen_desc,
            'POLYPHEN_SCORE': polyphen_score,
            'CADD_PHRED': anotacion.get('CADD_phred'),
            'CADD_RAW': anotacion.get('CADD_raw'),
            'NUM_SAMPLES_WITH_MUTATION': len(samples_with_mutation),
            'MEAN_VAF': statistics.mean([i[1] for i in samples_with_mutation]),
            'SAMPLES_WITH_MUTATION': ';'.join([i[0] for i in samples_with_mutation]),
            'COSMIC_MATCH': ';'.join([i[0] for i in cosmic_match]),
            'HEURISTIC': 'Y' if heuristic_filter(cosmic_match) else 'N',
            'PREVIOUSLY_IDENTIFIED': 'Y' if is_previously_identified() else 'N',
            'WHITELIST': 'Y' if is_whitelisted else 'N',
            'INTERNALLY_IDENTIFIED': internally_identified.get(var_id),
            'ARTIFACT': 'Y' if is_artifact() else 'N',
            'VAF_2': 'Y' if vaf >= 0.02 else 'N',
            'VAF_1': 'Y' if vaf >= 0.01 else 'N',
            'VAF_0': 'Y' if vaf < 0.01 else 'N',
            'VAF_10': 'Y' if vaf >= 0.1 else 'N',
            'VARSOME': create_varsome_link()
        }
        return res

    def print_item(item):
        nonlocal header_flag
        if header_flag:
            for field in item.keys():
                print(field, end='\t')
            header_flag = False
            print()

        for field in item.values():
            print(field, end='\t')
        print()

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

    def create_varsome_link():
        return f"https://varsome.com/variant/hg38/{vcf_chrom[3:]}-{vcf_pos}-{vcf_ref}-{vcf_alt}?annotation-mode=somatic&tissue-type=blod)"

    bypass_cosmic = False  # argumento para saltarse la comprobacion de COSMIC. SOLO DEBUG

    try:
        VAF_THRESHOLD = config.getfloat('PARSING', 'VAF_THRESHOLD')
        MAX_MAF_THRESHOLD = config.getfloat('PARSING', 'MAX_MAF_THRESHOLD')
        MAX_NUM_SAMPLES_GERMINAL_VAF = config.getint('PARSING', 'MAX_NUM_SAMPLES_GERMINAL_VAF')
        VAF_THRESHOLD_EXCEPTION_GENES = config.get('PARSING', 'VAF_THRESHOLD_EXCEPTION_GENES').split(',')

    except NoOptionError as e:
        print(f'ERROR: {e}')

    db = chip_cfg.connect_db(config)
    cosmic_table = get_cosmic() if not bypass_cosmic else {}
    whitelist_ncl = get_whitelist_ncl()
    whitelist_aa = get_whitelist_aa()
    previously_identified = get_previously_identified()
    internally_identified = get_internally_identified()
    # canonical_transcripts = get_canonical_transcripts()
    artifacts = get_artifacts()

    db_connection = TargetPanelMySQLConnector(config)
    panel = TargetPanel(config, db_connection)
    panel_genes = panel.genes()
    canonical_transcripts = panel.transcripts()

    header_flag = True

    path_outfile = os.path.split(tsv_outfile)[0]
    os.makedirs(path_outfile, exist_ok=True)
    with open(tsv_outfile, 'wt') as out_file:
        driver_writer = csv.writer(out_file, delimiter='\t')
        with gzip.open(annotated_vcf, 'rt') as vcf_file:
            for line_chunk in vcf_file:
                if line_chunk.startswith('##INFO=<ID=CSQ'):
                    # Extract CSQ header from vcf
                    csq_header = extract_CSQ_header(line_chunk)
                elif line_chunk.startswith('#CHROM'):
                    sample_names = get_sample_name(line_chunk.split()[9:], translation_filename)
                    num_samples = len(sample_names)
                    MAX_NUM_SAMPLES_WITH_VAR = math.ceil(
                        num_samples * config.getfloat('PARSING', 'TOTAL_SAMPLES_FRACTION'))
                elif not line_chunk.startswith('#'):
                    vcf_chrom, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, *vcf_samples \
                        = line_chunk.strip().split('\t')

                    variante = get_mutation(vcf_alt)
                    var_id = f"""{vcf_chrom}_{vcf_pos}_{vcf_ref}_{variante}"""

                    anotacion = get_annotation()

                    if anotacion and is_deleterious(anotacion):
                        samples_gt = get_sample_data([vcf_ref, *vcf_alt.split(',')], sample_names, vcf_format,
                                                     vcf_samples)

                        samples_with_mutation = get_samples_with_mutation(samples_gt)
                        samples_with_vaf_like_germinal = get_samples_with_vaf_like_germinal(samples_with_mutation)

                        is_whitelisted = is_whitelist_ncl(var_id) or is_whitelist_aa()
                        if is_whitelisted or seems_somatic():
                            for sample_id, vaf in samples_with_mutation:
                                ok_vaf = (vaf >= VAF_THRESHOLD or is_in_vaf_threshold_exception_gene()) \
                                    if not bypass_vaf_filter else True
                                if not is_vaf_germ(vaf, nocall_as_germ=True) and ok_vaf:
                                    item = create_printing_item()
                                    if header_flag:
                                        driver_writer.writerow(item.keys())
                                        header_flag = False

                                    driver_writer.writerow(item.values())
    db.close()
