import os
from gold.gsuite import GSuiteConstants
from gold.gsuite.GSuite import GSuite
from gold.gsuite.GSuiteComposer import composeToFile
from gold.gsuite.GSuiteFunctions import getTitleWithSuffixReplaced, writeGSuiteHiddenTrackStorageHtml
from gold.gsuite.GSuiteTrack import GalaxyGSuiteTrack, GSuiteTrack
from quick.extra.ProgressViewer import ProgressViewer
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.clustering.RsidMapper import RsidMapper


class GSuiteSumstatConversionTool(GeneralGuiTool):

    GSUITE_ALLOWED_FILE_FORMATS = [GSuiteConstants.PRIMARY, GSuiteConstants.UNKNOWN]
    GSUITE_ALLOWED_LOCATIONS = [GSuiteConstants.LOCAL]

    HISTORY_PROGRESS_TITLE = 'Progress'
    HISTORY_HIDDEN_TRACK_STORAGE = 'GSuite track storage'
    HISTORY_ERROR_TITLE = 'GSuite - files that failed manipulation'

    HG38 = 'hg38'
    HG19 = 'hg19'

    @staticmethod
    def getToolName():
        return "Convert sumstat files in GSuite to valued point tracks"

    @staticmethod
    def getInputBoxNames():
        return [
            ('Select GSuite from history', 'gSuite'),
            ('Filter rows by value', 'hasFilter'),
            ('Select effect size threshold to filter by', 'valFilter'),
            ('Log-transform values: -log(val)', 'logTransform'),
            ('Select reference genome for uniform mapping', 'refGenome')
        ]

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxHasFilter(choices):
        return False

    @staticmethod
    def getOptionsBoxValFilter(choices):
        if choices.hasFilter:
            return '0.0001'

    @staticmethod
    def getOptionsBoxLogTransform(choices):
        return False

    @staticmethod
    def getOptionsBoxRefGenome(choices):
        return [
            GSuiteSumstatConversionTool.HG38,
            GSuiteSumstatConversionTool.HG19
        ]

    @staticmethod
    def getInfoForOptionsBoxLogTransform(choices):
        return 'This transformation is done for p-values of SNPs in tracks extracted from the GWAS Catalog. Depending' \
               ' on the analysis the tracks will be used for, it can be an advantage to convert the values. ' \
               'One example is if the p-values are to be used in a comparison with Pearson correlation coefficient.'

    @classmethod
    def getExtraHistElements(cls, choices):
        from quick.webtools.GeneralGuiTool import HistElement
        from gold.gsuite.GSuiteConstants import GSUITE_SUFFIX, GSUITE_STORAGE_SUFFIX

        fileList = [HistElement(cls.HISTORY_ERROR_TITLE, GSUITE_SUFFIX)]
        fileList += [HistElement(cls.HISTORY_PROGRESS_TITLE, 'customhtml')]

        if choices.gSuite:
            try:
                fileList.append(HistElement(cls.HISTORY_HIDDEN_TRACK_STORAGE, GSUITE_STORAGE_SUFFIX, hidden=True))
                return fileList
            except:
                pass

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):
        from quick.webtools.clustering.GSuitePrimaryTrackModifier import GSuitePrimaryTrackModifier

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        outGSuite = GSuite()
        errorGSuite = GSuite()

        progressViewer = ProgressViewer([('Manipulate tracks', gSuite.numTracks() + 24)],
                                        cls.extraGalaxyFn[cls.HISTORY_PROGRESS_TITLE])

        rsidMap = RsidMapper.createRsidMappingFromStaticFiles(progressViewer, choices.refGenome)

        hiddenStorageFn = cls.extraGalaxyFn[cls.HISTORY_HIDDEN_TRACK_STORAGE]

        for track in gSuite.allTracks():
            fileName = cls.getFilenameWithGTrackSuffix(track.path)
            title = getTitleWithSuffixReplaced(track.title, 'gtrack')

            try:
                uri = GalaxyGSuiteTrack.generateURI(
                    galaxyFn=hiddenStorageFn,
                    extraFileName=fileName,
                    suffix='gtrack'
                )

                gSuiteTrack = GSuiteTrack(
                    uri,
                    title=title,
                    genome=track.genome,
                    trackType='valued points',
                    attributes=track.attributes
                )

                trackFn = gSuiteTrack.path
                if choices.hasFilter:
                    GSuitePrimaryTrackModifier.convertSumstat(track.path, trackFn, rsidMap,
                                                              choices.logTransform, float(choices.valFilter))
                else:
                    GSuitePrimaryTrackModifier.convertSumstat(track.path, trackFn, rsidMap, choices.logTransform)

                outGSuite.addTrack(gSuiteTrack)

            except Exception as e:
                track.comment = 'An error occurred for the following track: ' + str(e)
                errorGSuite.addTrack(track)

            progressViewer.update()

        outGSuite.setGenomeOfAllTracks(choices.refGenome)
        composeToFile(outGSuite, galaxyFn)
        composeToFile(errorGSuite, cls.extraGalaxyFn[cls.HISTORY_ERROR_TITLE])
        writeGSuiteHiddenTrackStorageHtml(hiddenStorageFn)

    @classmethod
    def validateAndReturnErrors(cls, choices=None):

        errorStr = GeneralGuiTool._checkGSuiteFile(choices.gSuite)
        if errorStr:
            return errorStr

        if choices.hasFilter:
            try:
                float(choices.valFilter)
            except:
                return 'Filter threshold is not a number'

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

    @staticmethod
    def getOutputFormat(choices=None):
        return 'gsuite'

    @classmethod
    def getFilenameWithGTrackSuffix(cls, path):
        filename = os.path.basename(path)
        prefix, suffix = os.path.splitext(filename)

        if suffix:
            return prefix + '.gtrack'
        else:
            return filename + '.gtrack'

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header("Convert sumstat primary tracks")
        core.paragraph("Tool for conversion of GSuite with primary tracks that have the sumstat format [1], to a "
                       "GSuite of valued point tracks. The resulting valued point tracks contain columns for "
                       "chromosome number (seqid), chromosome position (start), SNP rsid (id) and P or Z values "
                       "(value). If the sumstat contains Z-values, no conversion of the value is needed. However, if "
                       "it has P-values, it is recommended to log transform them. A threshold for the value can be set,"
                       " and all track elements with values <= the given threshold will be kept in the new GSuite.")
        core.divider()
        core.smallHeader("Sumstat primary track requirements")
        core.paragraph("The tool assume that all primary tracks are tab separated. In addition, it requires the "
                       "following column headers in the original tracks to be present: snp and either p or z."
                       "The 'snp' column contains the rsids of all SNPs in the sumstat files, while the 'p' or 'z' "
                       "column contains the effect of the SNP. "

                       "Note that if the rsid of a track element is not present in the files within the "
                       "HyperBrowser that map rsids to positions on the given reference genome, the track element and "
                       "row will be dropped.")
        core.divider()
        core.smallHeader("Example use")
        core.paragraph("The sumstat files, as generated from the software of [2], does not contain positional "
                       "information of the summary statistic SNPs. This tool can be used to convert a collection "
                       "(GSuite) of sumstat files into a GSuite of valued point tracks. This GSuite can further "
                       "be used in HyperBrowser analyses.")
        core.paragraph("[1]: https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format<br>"
                       "[2]: Bulik-Sullivan et al. (2015). An atlas of genetic correlations across human diseases and "
                       "traits. In: Nature Genetics 47.11, pp.1236-1241.")
        core.divider()
        core.paragraph("<b>NB:</b> Tool takes some time to run, as it reads in a full rsid-position mapping for all dbSNP rsids "
                       "in the chosen reference genome.")
        return str(core)
