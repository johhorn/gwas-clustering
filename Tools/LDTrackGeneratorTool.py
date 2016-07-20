from gold.application.HBAPI import doAnalysis
from gold.gsuite import GSuiteConstants
from gold.track.Track import Track
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.mixin.DebugMixin import DebugMixin


class LDTrackGeneratorTool(GeneralGuiTool, DebugMixin):

    GSUITE_ALLOWED_FILE_FORMATS = [GSuiteConstants.PREPROCESSED]
    GSUITE_ALLOWED_LOCATIONS = [GSuiteConstants.LOCAL]

    DEFAULT_SELECT = '--- Select ---'

    @staticmethod
    def getToolName():
        return "Create a master LD graph track for a GSuite (hg38)"

    @staticmethod
    def getInputBoxNames():
        return [
            ('Select GSuite to generate LD track for', 'gSuite'),
            ('Edges are undirected', 'isUndirected'),
            ('Select r<sup>2</sup> treshold', 'rsquare')
        ] + DebugMixin.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxIsUndirected(choices):
        return True

    @staticmethod
    def getOptionsBoxRsquare(choices):
        return [LDTrackGeneratorTool.DEFAULT_SELECT, '1.0', '0.9', '0.8', '0.7', '0']

    @staticmethod
    def isDebugMode():
        return True

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):
        from gold.description.AnalysisDefHandler import AnalysisSpec
        from quick.statistic.UniquePointTrackStat import UniquePointTrackStat
        from quick.application.UserBinSource import GlobalBinSource

        cls._setDebugModeIfSelected(choices)

        # Analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)
        analysisSpec = AnalysisSpec(UniquePointTrackStat)

        # Get unique rsids for whole GSuite
        rsids = set()
        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            result = doAnalysis(analysisSpec, analysisBins, [track]).getGlobalResult()
            if 'Result' in result:
                rsids.update(result['Result'])

        # Create linked point track
        cls.createLinkedPointTrack(rsids, str(choices.isUndirected), galaxyFn, float(choices.rsquare))

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
            allowedFileFormats=cls.GSUITE_ALLOWED_FILE_FORMATS,
            allowedLocations=cls.GSUITE_ALLOWED_LOCATIONS
        )

        if errorStr:
            return errorStr

        if choices.rsquare == LDTrackGeneratorTool.DEFAULT_SELECT:
            return 'Please select a threshold of r<sup>2<sup>'

    @staticmethod
    def getOutputFormat(choices=None):
        return 'gtrack'

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def createLinkedPointTrack(cls, rsids, isUndirected, trackFn, r2):
        from quick.webtools.clustering.CreateLDTrack import CreateLDTrack
        from quick.application.ExternalTrackManager import ExternalTrackManager

        # Create file for GTrack
        galaxyTN = ExternalTrackManager.constructGalaxyTnFromSuitedFn(trackFn, fileEnding='gtrack', name='ld_graph')
        fn = ExternalTrackManager.extractFnFromGalaxyTN(galaxyTN)
        f = open(fn, 'w')

        # Get LD information and create linked point track
        ldDict = CreateLDTrack.getLDDict(r2)
        expansionDict = CreateLDTrack.getExpansionDict(rsids, ldDict)
        f.write(CreateLDTrack.formatLinkedPointTrack(expansionDict, isUndirected))

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header("Create master LD track from GSuite")
        core.paragraph("This tool takes in a GSuite, and generate a master track of LD information, across all the "
                       "original primary tracks. The new LD track is formatted as a linked point track. This track "
                       "represents an LD graph where the SNPs are nodes, and variants in LD are connected through "
                       "edges weighted by r<sup>2</sup>.")
        core.paragraph("If 'Edges are undirected' is checked, the graph edges will be undirected. The r<sup>2</sup> "
                       "threshold sets the lower limit of LD correlation. r<sup>2</sup> >= 0.8 is the consensus "
                       "threshold for a significant LD correlation. Only variants with r<sup>2</sup> >= the chosen "
                       "threshold will be included in the linked point track.")
        core.divider()
        core.smallHeader("Primary track requirements")
        core.paragraph("The primary tracks need a column header of 'snps'. This column must contain rsids, which "
                       "represent the different track elements.")
        core.divider()
        core.smallHeader("Limitations")
        core.paragraph("The tool makes use of LD information from a static file within the HyperBrowser to expand the "
                       "primary tracks with LD information. This file was generated with all significant SNPs from the "
                       "GWAS Catalog (18660 from last extraction). If the GSuite tracks contain rsids that are not "
                       "present in the LD master file, the track element will not be part of the new linked point "
                       "tracks.")
        core.paragraph("The master LD file was created for the CEU population. Thus, this tool is not appropriate if "
                       "information of LD between SNPs of other populations is wanted.")
        core.paragraph("The resulting GSuite is mapped to hg38.")
        core.divider()
        core.smallHeader("Example use")
        core.paragraph("Can be used to create a master linked point track for a GWAS Catalog GSuite with LD variants. "
                       "This track can further be utilized to explore the LD structure through the tool "
                       "'Empirical exploration of LD tracks'. Alternatively, one can use it with the original GSuite "
                       "in the tools 'Direct match with LD expansion' and 'Bipartite matching of valued point tracks'. ")
        return str(core)
