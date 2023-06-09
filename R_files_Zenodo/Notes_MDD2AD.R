
select distinct subject_name, subject_cui, count(*) as cnt from eidos where predicate = 'Influence' and object_cui in ('C0002395', 'C0494463', 'C0276496') group by subject_name, subject_cui order by subject_name asc;

select distinct subject_name, subject_cui, count(*) as cnt from eidos where predicate = 'Influence' and object_cui in ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')  group by subject_name, subject_cui order by subject_name asc;


kgConfounders <- c('Amyloidosis', 'Amyotrophic_Lateral_Sclerosis', 'Anemia', 'Apraxias', 'Atherosclerosis', 'Atrial Fibrilation', 'Brain hemorrhage', 'Cerebral_Amyloid_Angiopathy', 'Cerebral_atrophy', 'Cerebrovascular_accident', 'Chronic_infectious_disease', 'Chronic_Obstructive_Airway_Disease', 'Congestive heart failure', 'Deglutition Disorders', 'Diabetes_Mellitus', 'Diabetes_Mellitus,_Non-insulin-Dependent', 'Encephalitis', 'Encephalopathies', 'Falls', 'Frontotemporal_dementia', 'Hypercholeresterolemia', 'Hyperglycemia', 'Hyperinsulinism', 'Hypertensive_disease', 'Hypoglycemia', 'Hypotension', 'Hypotension,_orthostatic', 'Increase in blood pressure', 'INFECTION RECURRENT', 'Immune_response', 'Inflammatory_response', 'Insulin_Resistance', 'Ischemic_stroke', 'Leukoencephalopathy', 'Low tension glaucoma', 'Malnutrition', 'Migraine_Disorders', 'Myocardial infarction', 'Nerve_Degeneration', 'Neurofibrillary_degeneration_(morphologic_abnormality)', 'Non-alcoholic_Fatty_Liver_Disease', 'Obesity', 'Osteoporosis', 'Overweight', 'Parkinsonian_Disorders', 'Periodontal_Diseases', 'Pneumonia', 'Psychotic_Disorders', 'Rickets', 'Senile_Plaques', 'Sleep_Apnea_Syndrome', 'Sleep_Apnea,_Obstructive', 'Subarachnoid Hemorrhage', 'Tauopathies', 'Unconscious State', 'Vitamin D Deficiency')

confounders[,4]

intersect(tolower(kgConfounders), confounders[,4])


length(kgConfounders)



###### bacs is a very strange category!

### exclusions!
# 'medd' # placebos, drug delivery systems, probs, self-help devices
# 'genf' # cellular processes - splicing
# 'imft' # vaccines, antibodies, epitopes
# 'chvs' # derivatives, lead compound, monomer, particle, dimer, macromolecule
# 'orga' # sex, ability, female, gender, genetic predisposition to disease, adiposity
# 'clna' 
# 'bdsu' # cerebrospinal fluid, plasma, platelet rich plasma
# 'nnon' # RNA (untranslated), transcript, RNA (circular)
# 'nusq' # SNP, transcriptome, DNA sequence
# 'chem' # ligands, acids, chemicals
# 'diap' # MRI, electroencephalography, PET, neuroimaging
# 'chvf' # agent, inhibitors, antagonists, molecular target
# 'sbst' # constituents, players, substance, liquid substance, breath, catabolites
# 'bodm' # injection, capsule, dendrimes
# 'resd' # study models, biological models
# drdd' # skin patch
###
# ('medd', 'genf', 'imft', 'chvs', 'orga', 'clna', 'bdsu', 'nnon', 'nusq', 'chem', 'diap', 'chvf', 'sbst', 'bodm', 'resd', 'drdd')


# drug semtypes
# 'phsu' # cholinesterase inhibitors, donepezil, acatylcholinesterase inhibitors, memantine
# 'orch'
# 'hops' # lead, mercury, tobacco, DDT
# 'food'
# 'irda' # curcumin, acrolein, indicators
# 'clnd' # donepezil
# 'vita' # daizdein, lipoate, vitamin d, thiamine monophosphate
# 'eico' # arachidonic acid, dinoprostone
# 'opco' # phosphorylcholine, diphosphonates, phosphoric acid esters
# #########
#
# ('phsu', 'orch', 'hops', 'food', 'irda', 'clnd', 'vita', 'eico', 'opco')

# genes/enzymes
# 'aapp'
# 'gngm'
# 'horm'
# 'bacs'
# 'lipd'
# 'elii' # aluminum, ions, mercury, carbon
# 'inch'
# 'antb' # stretozocin, sirolimus, antibiotics
# 'carb '# galactose, ascorbic acid, polysacchardies, blood glucose, carbohydrates
# 'nsba' # neurotransmitters, also dopamine, amphetamine, serotonin
#######
# ('aapp', 'gngm', 'horm', 'bacs', 'lipd', 'elii', 'inch', 'antb', 'carb', 'nsba')


####################
####################
# procedures, interventions
# 'topp'
# 'dora' # <- walking, exercise, reading, walking
####################
####################
### habits
# 'inbe' # smoking, alcohol consumption, habits, participant, defensiveness, desires, dietary habits
# ##################
####################
### diseases/syndromes
# 'dsyn'
# 'fndg' # impaired cognition, stress, high fat diet
# 'patf'
# 'sosy'
# 'mobd' # dementia, memory impairment, schizophrenia, ...
# 'strd'
# 'inpo' # toxic effect, TBI (!), brain injury, wounds and injuries
# 'lbtr' # <- bone density
# 'virs' # herpesvirus 1 (human), herpesviridae, 
# 'neop' # malignant neoplasms, malignant fibrous histiocytoma
# 'bact' # 
# 'invt' # parasites, tetameres, toxoplasma gondii
# 'orgm' # infectious agent, pathogenic organism
# 'cgab' # down syndrome, congenital abnormality
# 'rich'
# 'topp'
# 'dora'
# 'acab' # neurofibrillary tangles
# 'eehu' # air pollution
####
# ('dsyn', 'fndg', 'patf, 'sosy', 'mobd', 'strd', 'inpo', 'lbtr', 'virs', 'neop', 'bact', 'invt', 'orgm', 'cgab', 'rich', 'topp', 'dora', 'acab', 'inbe'', 'eehu')




# with xy as (select distinct subject_semtype, subject_name from causalpredications where predicate IN ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS') AND object_cui IN ('C0002395', 'C0494463', 'C0276496')),
#xx as (select distinct subject_semtype, subject_name from causalpredications where predicate IN ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS') AND object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517'))
#select distinct subject_semtype, count(*) cnt from xy group by subject_semtype order by cnt DESC;

# subjects by count 
# select distinct subject_name, count(*) cnt from causalpredications where predicate IN ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS', 'INHIBITS', 'STIMULATES', 'AFFECTS', 'DISRUPTS') and object_cui IN ('C0002395', 'C0494463', 'C0276496') and subject_semtype = 'rich' group by subject_name order by cnt desc;


# ('C0002395', 'C0494463', 'C0276496')


with xy as (select distinct subject_semtype, subject_name from causalpredications where predicate IN ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS') AND object_cui IN ('C0002395', 'C0494463', 'C0276496')),
xx as (select distinct subject_semtype, subject_name from causalpredications where predicate IN ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS') AND object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517'))
select distinct subject_semtype, count(*) cnt from xy group by subject_semtype order by cnt DESC;



