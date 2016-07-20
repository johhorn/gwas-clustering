from quick.statistic.UniquePointTrackStat import UniquePointTrackStat
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from gold.application.HBAPI import doAnalysis, AnalysisSpec, GlobalBinSource
from gold.application.HBAPI import Track
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin
from gold.result.HtmlCore import HtmlCore


class PrintingTool(GeneralGuiTool, DebugMixin):

    # Printing cases related to external LD expansion
    DEFAULT = '--- Select ---'
    PRINT_SNPS = 'Print all unique SNPs in GSuite'
    PRINT_DISEASE = 'Print unique SNPs for each track in GSuite'

    @staticmethod
    def getToolName():
        return "Print SNPs of GSuite"

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def getInputBoxNames(cls):
        return [
                    ('Select GSuite from history', 'gSuite'),
                    ('Select empirical case', 'getEmpiricalCase')
               ] + cls.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxGSuite():
            return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxGetEmpiricalCase(choices):
        return [
            PrintingTool.DEFAULT,
            PrintingTool.PRINT_SNPS,
            PrintingTool.PRINT_DISEASE
        ]

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):

        cls._setDebugModeIfSelected(choices)

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)

        core = HtmlCore()
        core.divBegin(style=CommonClusteringFunctions.HTML_STYLE)

        if choices.getEmpiricalCase == PrintingTool.PRINT_SNPS:
            core.header(PrintingTool.PRINT_SNPS)
            core.divider(True)
            cls.printUniqueRsIDs(analysisBins, gSuite, core)
        elif choices.getEmpiricalCase == PrintingTool.PRINT_DISEASE:
            core.header(PrintingTool.PRINT_DISEASE)
            core.divider(True)
            cls.printDiseaseRsIDs(analysisBins, gSuite, core)

        core.divEnd()
        print core

    @staticmethod
    def isDebugMode():
        return True

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @classmethod
    def printUniqueRsIDs(cls, analysisBins, gSuite, core):
        rsIDs = set()
        analysisSpec = AnalysisSpec(UniquePointTrackStat)

        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            result = doAnalysis(analysisSpec, analysisBins, [track])
            if 'Result' in result.getGlobalResult():
                observations = result.getGlobalResult()['Result']
                rsIDs.update(observations)

        core.smallHeader('Printing snps from GSuite of ' + str(gSuite.numTracks()) + ' tracks')
        for rsid in rsIDs:
            core.line(rsid)

    @classmethod
    def printDiseaseRsIDs(cls, analysisBins, gSuite, core):
        analysisSpec = AnalysisSpec(UniquePointTrackStat)

        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            result = doAnalysis(analysisSpec, analysisBins, [track])

            core.smallHeader('<br><br><br>' + gSuiteTrack.title + '<br><br>')
            if 'Result' in result.getGlobalResult():
                observations = result.getGlobalResult()['Result']
                for snp in observations:
                    core.line(snp)

    @classmethod
    def validateAndReturnErrors(cls, choices=None):

        errorStr = GeneralGuiTool._checkGSuiteFile(choices.gSuite)
        if errorStr:
            return errorStr

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        if gSuite.numTracks() == 0:
            return 'Please select a GSuite file with at least one track'

        errorStr = cls._checkGSuiteRequirements(
            gSuite,
            allowedFileFormats=CommonClusteringFunctions.GSUITE_ALLOWED_FILE_FORMATS,
            allowedLocations=CommonClusteringFunctions.GSUITE_ALLOWED_LOCATIONS
        )

        if errorStr:
            return errorStr

        if choices.getEmpiricalCase == cls.DEFAULT:
            return 'Please select a printing option'

    @staticmethod
    def getToolDescription():
        core = HtmlCore()
        core.header('Print SNPs in valued point tracks')
        core.paragraph('This tool prints out the unique SNP rsids for all track elements in a GSuite.'
                       'It has been made for point tracks with a column header of "snps", whose rows are SNP rsids. '
                       'This is typically a track extracted from the GWAS Catalog within the HyperBrowser.')
        core.divider()
        core.smallHeader('Example usage')
        core.paragraph('This tool is useful when a GSuite of SNPs need to be externally processed. '
                       'For instance, if one wish to find variants in LD from a given GSuite, using the Ensembl '
                       'Variation API. One can print out all unique rsids of the GSuite and extract variant '
                       'information, which in turn can be uploaded as a new track.')
        return str(core)
