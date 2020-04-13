                amplimer.hit_info = amplimer.hit_info.replace('\t', '').strip().split('\n')
                amplimer.hit_info = [list(filter(None, i.split(' '))) for i in amplimer.hit_info]

                amplicon = {'hit_no': index}
                amplicon['template_name'] = amplimer.hit_info[0][0]
                amplicon['full_template_seq'] = self.template#seq[seq_idx[amplicon['template_name']]]
                amplicon['forward_oligo'] = {'primer_seq': Seq(self._collapsed_iupac(amplimer.hit_info[1][0]), alphabet = IUPAC.ambiguous_dna),
                                        'strand': amplimer.hit_info[1][2],
                                        'position': int(amplimer.hit_info[1][5]) - 1,
                                        'mismatches': amplimer.hit_info[1][7]}
                amplicon['reverse_oligo'] = {'primer_seq': Seq(self._collapsed_iupac(amplimer.hit_info[2][0]), alphabet = IUPAC.ambiguous_dna).reverse_complement(),
                                        'strand': amplimer.hit_info[2][2],
                                        'position': len(amplicon['full_template_seq'].seq) - int(amplimer.hit_info[2][5].replace('[', '').replace(']', '')) + 1,
                                        'mismatches': int(amplimer.hit_info[2][7])}
                amplicon['target_template_seq'] = amplicon['full_template_seq'][amplicon['forward_oligo']['position']:amplicon['reverse_oligo']['position']]
                amplicon['target_template_seq'].description = f"| {primerpair}~~~Amplimer {index}~~~{amplicon['forward_oligo']['position']}:{amplicon['reverse_oligo']['position']}"
                amplicon['padded_primer_forward'] = SeqRecord(Seq(f"{amplicon['forward_oligo']['primer_seq']}{'-'*(len(amplicon['target_template_seq'])-len(amplicon['forward_oligo']['primer_seq']))}", alphabet = IUPAC.ambiguous_dna),#-len(amplimer['reverse_oligo']['primer_seq']))}{amplimer['reverse_oligo']['primer_seq']}"
                                        id = f'forward_oligo_{str(primerpair)}',
                                        description = str(primerpair))
                amplicon['padded_primer_reverse'] = SeqRecord(Seq(f"{'-'*(len(amplicon['target_template_seq'])-len(amplicon['reverse_oligo']['primer_seq']))}{amplicon['reverse_oligo']['primer_seq']}", alphabet = IUPAC.ambiguous_dna),#-len(amplimer['reverse_oligo']['primer_seq']))}{amplimer['reverse_oligo']['primer_seq']}"
                                        id = f'reverse_oligo_{str(primerpair)}',
                                        description = str(primerpair))
                results_dfs_list.append(amplicon, index=[primerpair_name])
