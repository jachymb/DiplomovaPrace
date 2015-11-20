import pickle
import gzip
import xml.etree.ElementTree as ElementTree
import re
import operator
from dp.utils import debug, parseFasta
from pathlib import Path
from collections import defaultdict
from urllib.request import urlopen, Request

__all__ = ["GeneFactory", "Gene"]

PICKLEDIR = Path('data/serialized_genes')
XMLDIR = Path('data/xml')
SS = Path('data/ss.txt')

MAPPING = {
        'A' : 'ala',
        'R' : 'arg',
        'N' : 'asn',
        'D' : 'asp',
        'B' : 'asx', # asparagine or aspartic acid
        'C' : 'cys',
        'E' : 'glu',
        'Q' : 'gln',
        'Z' : 'glx',
        'G' : 'gly',
        'H' : 'his',
        'I' : 'ile',
        'J' : 'xle', #Leucine or Isoleucine
        'L' : 'leu',
        'K' : 'lys',
        'M' : 'met',
        'O' : 'pyl', # Pyrrolysine
        'F' : 'phe',
        'P' : 'pro',
        'S' : 'ser',
        'T' : 'thr',
        'W' : 'trp',
        'X' :  None, # Unknown
        'Y' : 'tyr',
        'U' : 'sec', # Selenocysteine
        'V' : 'val',
        'Z' : 'glx', # Glutamine or glutamic acid
        }

MAPPING_SS = {
        'T': 'turn',
        'G': 'helix_3_10',
        'S': 'bend',
        'H': 'alpha_helix',
        'E': 'beta_stand',
        'B': 'beta_bridge',
        'I': 'pi_helix',
        ' ': None}

IGNORE = (None, 'sec', 'pyl')
#IGNORE = (None,)

class GeneFactory:
    openGenes = []

    def __init__(self, deserialize = True):
        self.deserialize = deserialize

    def getGene(self, fullName):
        fullName = Gene.canonicalName(fullName)
        name = Gene.canonicalName(fullName, False)
        serializedFileName = Gene.serializedFileName(fullName)
        xmlFileName = Gene.xmlname(name)
        while name in self.openGenes:
            time.sleep(0.1)
        if self.deserialize and serializedFileName.is_file() and serializedFileName.stat().st_mtime > xmlFileName.stat().st_mtime:
            #debug("Deserializing data for gene " + fullName)
            with serializedFileName.open('rb') as f:
                try:
                    data = pickle.load(f)
                    return Gene(*data)
                except EOFError:
                    pass

        self.openGenes.append(name)
        found = None
        for gene in Gene.fromXML(fullName):
            if gene.name == fullName: found = gene
        self.openGenes.remove(name)
        
        if found is None:
            raise KeyError(fullName)
        return found

class Gene:
    """ This class represents a single gene and the protein it encodes. """

    DISTANCE_RESOLUTION = 2
    MAX_DISTANCE = 14

    XML_URL = "http://www.rcsb.org/pdb/files/%s.xml.gz"
    recalculateDists = False

    ss = dict(parseFasta(SS))
    __slots__ = ["name", "structure", "sequence", "distances", "secstr"]
    def __init__(self, name, sequence, structure, distances = None, secstr = None):
        self.name = name
        self.structure = structure
        self.sequence = sequence
        calcDists = (distances is None or self.recalculateDists)
        self.distances = tuple(self.getDistances() if calcDists else distances)
        if calcDists:
            self.dump()
        if secstr is None:
            self.getSecStr()
        else:
            self.secstr = secstr
        #debug('Initialized gene '+name)

    def getSecStr(self):
        name, strand = self.name.split("_")
        name = name.upper()
        try:
            s = tuple([MAPPING[a] for a in self.ss["%s:%s:sequence" % (name, strand)]])
            if s != self.sequence:
                debug("WARNING: Different sequences for %s" % (self.name,))
                debug(str(s))
                debug(str(self.sequence))
            self.secstr = self.ss["%s:%s:secstr" % (name, strand)]
        except KeyError:
            debug("WARNING: missing secondary structure info %s" % (self.name,))
            self.secstr = None
        self.dump()


    @staticmethod
    def canonicalName(name, withStrand = True):
        name, strand = name.split("_")
        name = name.lower()
        if withStrand:
            return name+"_"+strand
        else:
            return name

    @staticmethod
    def xmlname(name):
        name = name.lower()
        return XMLDIR / name[1:3] / (name + ".xml.gz")

    @classmethod
    def fromXML(cls, fullName):
        """Generates protein data from stored XML file. Downloads it if not present."""
        name = Gene.canonicalName(fullName, False)
        fname = cls.xmlname(name)
        sequence = []
        structure = defaultdict(dict) # structure[seq_id][seq_position] == coordinates
        sequences = defaultdict(dict)
        strand2seq = {}

        # Download data if not stored on disc
        if not fname.is_file():
            debug("Downloading file for gene %s... " % name, False)
            req = Request(cls.XML_URL % name, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.47 Safari/537.36'})
            with fname.open('wb') as f:
                f.write(urlopen(req).read())
            #debug("Done.")

        # Parse the file
        debug("Parsing XML for gene %s... " % name)

        tag_stack = []
        elem_stack = []

        pp_seq_id = None
        with gzip.open(str(fname), 'r') as f:
            #fcntl.lockf(f.fileno(), fcntl.LOCK_EX)
            doc = ElementTree.iterparse(f, ('start','end'))

            #next(doc) # Skip root element
            _, root = next(doc)
            XMLNS = root.tag[:root.tag.find("}")+1]

            path_atom = [XMLNS + 'atom_siteCategory', XMLNS + 'atom_site']
            path_seq = [XMLNS + 'entity_poly_seqCategory', XMLNS + 'entity_poly_seq']
            path_poly = [XMLNS + 'entity_polyCategory', XMLNS + 'entity_poly']

            whitespace = re.compile(r'\s+')
            for event, elem in doc:
                if event == 'start':
                    tag_stack.append(elem.tag)
                    elem_stack.append(elem)
                elif event == 'end':
                    if tag_stack == path_atom:
                        # elem[11] = <PDBx:label_atom_id>, elem[5] = <PDBx:auth_atom_id>
                        if elem[11].text == 'CA' and elem[5].text == 'CA':
                            # elem[13] = <PDBx:label_entity_id>
                            seq_id = elem[13].text

                            # elem[14] = <PDBx:label_seq_id>, elem[{1,2,3}] = <PDBx:Cartn_{xyz}>

                            label_seq_id = elem[14].text
                            if label_seq_id is not None:
                                coordinates = [float(elem[i].text) for i in (1,2,3)]
                                seq_pos = int(elem[14].text) - 1
                                structure[seq_id][seq_pos] = tuple(coordinates)
                        elem_stack[-2].remove(elem)

                    elif tag_stack == path_poly:
                        seq_type = elem.find(XMLNS + 'type').text
                        if seq_type.startswith('polypeptide'): # elem[5] = <PDBx:type>
                            seq = elem.find(XMLNS + 'pdbx_seq_one_letter_code_can').text
                            sequence = [MAPPING[c] for c in re.sub(whitespace, '', seq)]
                            seq_id = elem.attrib['entity_id']
                            strand_ids = elem.find(XMLNS + 'pdbx_strand_id').text
                            strand2seq.update({strand_id : seq_id for strand_id in strand_ids.split(',')})
                            sequences[seq_id] = sequence
                        elem_stack[-2].remove(elem)

                    #elif tag_stack == path_seq:
                        # entity_id="1" shoud be the polypeptide sequence, I hope, but I'm not sure. FIXME
                        # There can also be e.g. polynucleotide sequences
                    #    if elem.attrib['entity_id'] == '1':
                    #        sequence.append(elem.attrib['mon_id'].lower())
                    #    elem_stack[-2].remove(elem)


                    if tag_stack: tag_stack.pop()
                    if elem_stack: elem_stack.pop()

            #atomQuery = "./%satom_siteCategory/%satom_site[%slabel_atom_id='N'][%sauth_atom_id='N'][%slabel_entity_id='1']" % ((XMLNS,) * 5)
            #seqQuery = "./%sentity_poly_seqCategory/%sentity_poly_seq[@entity_id='1']" % ((XMLNS,) * 2)

        #debug("Done.")
        for strand, seq in strand2seq.items():
            yield cls(name+'_'+strand, tuple(sequences[seq]), tuple(structure[seq].items()))

    def getDistances(self):
        """Calculate distances among amino acid units. (Distance of N atom in amine.)"""
        for i, x in self.structure:
            for j, y in self.structure:
                if i <= j: continue
                dist = round(sum([x**2 for x in map(operator.sub, x, y)]) / self.DISTANCE_RESOLUTION) * self.DISTANCE_RESOLUTION
                if dist <= self.MAX_DISTANCE:
                    yield (i,j, dist)

    def resName(self, index):
        """Returns canonical name of residue at given position."""
        return "res_%d_%s" % (index+1, self.name.lower())

    def logicalRepresentation(self, postfix = None):
        """Generates the logical representation usable for learning."""

        results = [("proteinName", (self.name,))] # predicate head, parameters
        #results = []

        if postfix is None:
            postfix = self.name

        last = None
        # Residue sequence and sequence data
        for i, aminoAcid in enumerate(self.sequence):
            resName = self.resName(i)
            results.append(("res", (resName,)))
            if aminoAcid not in IGNORE:
                results.append(("residue", (resName, aminoAcid)))
            if self.secstr:
                name = MAPPING_SS[self.secstr[i]]
                if name is not None:
                    results.append(("secstr", (resName, name)))


            #if i > 0:
            #    results.append(("next", (self.resName(i-1), resName)))

        # Distance literals
        for i, j, dist in self.distances:
            #results.append(("dist", (self.resName(i), self.sequence[i], self.resName(j), self.sequence[j], dist)))
            results.append(("dist", (self.resName(i), self.resName(j), dist)))
                #results.sort()
        return ["%s(%s)" % (head, ", ".join(map(str, args))) for head, args in results] + self.backgroundKnowledge

    @staticmethod
    def serializedFileName(name):
        return PICKLEDIR / (name+".pickle")

    def dump(self):
        data = (self.name, self.sequence, self.structure, self.distances, self.secstr)
        with self.serializedFileName(self.name).open('wb') as f:
            pickle.dump(data, f)

