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


class GSuiteLiftOverTool(GeneralGuiTool):

    GSUITE_ALLOWED_FILE_FORMATS = [GSuiteConstants.PRIMARY, GSuiteConstants.UNKNOWN]
    GSUITE_ALLOWED_LOCATIONS = [GSuiteConstants.LOCAL]

    HISTORY_PROGRESS_TITLE = 'Progress'
    HISTORY_HIDDEN_TRACK_STORAGE = 'GSuite track storage'
    HISTORY_ERROR_TITLE = 'GSuite - files that failed manipulation'

    HG38 = 'hg38'
    HG19 = 'hg19'

    @staticmethod
    def getToolName():
        return "Genomic liftover for point track GSuites"

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def getInputBoxNames():
        return [
            ('Select GSuite from history', 'gSuite'),
            ('Select reference genome for uniform mapping', 'refGenome')
        ]

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxRefGenome(choices):
        return [
            GSuiteLiftOverTool.HG38,
            GSuiteLiftOverTool.HG19
        ]

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

        # Set analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        outGSuite = GSuite()
        errorGSuite = GSuite()
        progressViewer = ProgressViewer([('Manipulate tracks', gSuite.numTracks() + 24)],
                                        cls.extraGalaxyFn[cls.HISTORY_PROGRESS_TITLE])
        hiddenStorageFn = cls.extraGalaxyFn[cls.HISTORY_HIDDEN_TRACK_STORAGE]

        # Get rsID map for the chosen reference genome
        rsidMap = RsidMapper.createRsidMappingFromStaticFiles(progressViewer, choices.refGenome)

        # Lift over all tracks
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
                    trackType=track.trackType,
                    attributes=track.attributes
                )

                trackFn = gSuiteTrack.path
                GSuitePrimaryTrackModifier.liftOverGTrack(track.path, trackFn, rsidMap)
                outGSuite.addTrack(gSuiteTrack)

            except Exception as e:
                track.comment = 'An error occurred for the following track: ' + str(e)
                errorGSuite.addTrack(track)

            progressViewer.update()

        # Update reference genome of all tracks and write to file
        outGSuite.setGenomeOfAllTracks(choices.refGenome)
        composeToFile(outGSuite, galaxyFn)
        composeToFile(errorGSuite, cls.extraGalaxyFn[cls.HISTORY_ERROR_TITLE])
        writeGSuiteHiddenTrackStorageHtml(hiddenStorageFn)

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
        core.header("GSuite liftover")
        core.paragraph("Tool for lifting over all primary tracks within a GSuite to a specified reference genome. "
                       "Is made for point tracks (of SNPs) first and foremost, and depends on unique rsids to "
                       "identify each track element that should be lifted over.")
        core.divider()
        core.smallHeader("Primary track requirements")
        core.paragraph("The tool assume that all primary tracks are tab separated. In addition, it requires the "
                       "following column headers in the original tracks to be present: id, seqid, start. The 'id' "
                       "column contains rsid, the 'seqid' column contains chromosome number and 'start' column "
                       "contains chromosome position for all track elements. The columns of 'start' and 'seqid' will "
                       "be updated by the tool, with new positions relative to the chosen reference genome. All other "
                       "columns are left untouched. Note that if the rsid is not present in the files within the "
                       "HyperBrowser that map rsids to positions on the given reference genome, the track element and "
                       "row will be dropped.")
        core.divider()
        core.smallHeader("Example use")
        core.paragraph("If two GSuites of different reference genomes are to be compared, and positional information "
                       "is important for the analysis, one of the GSuites should be lifted over to the same reference "
                       "genome as the other.")
        core.paragraph("Another scenario where this tool is needed, is when primary tracks within a GSuite, or track "
                       "elements within primary tracks are mapped to different reference genomes, and positions should "
                       "be consistent across the entire GSuite. The tool will use the rsid column and lift over all "
                       "elements within the same GSuite to the chosen reference genome.")
        core.divider()
        core.paragraph("<b>NB:</b> Tool takes some time to run, as it reads in a full rsid-position mapping for all dbSNP rsids "
                       "in the chosen reference genome.")
        return str(core)
