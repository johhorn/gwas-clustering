import os
from gold.gsuite import GSuiteConstants
from gold.gsuite.GSuite import GSuite
from gold.gsuite.GSuiteComposer import composeToFile
from gold.gsuite.GSuiteFunctions import getTitleWithSuffixReplaced, writeGSuiteHiddenTrackStorageHtml
from gold.gsuite.GSuiteTrack import GalaxyGSuiteTrack, GSuiteTrack
from quick.extra.ProgressViewer import ProgressViewer
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.webtools.GeneralGuiTool import GeneralGuiTool


class LDGSuiteGeneratorTool(GeneralGuiTool):

    GSUITE_ALLOWED_FILE_FORMATS = [GSuiteConstants.PRIMARY, GSuiteConstants.UNKNOWN]
    GSUITE_ALLOWED_LOCATIONS = [GSuiteConstants.LOCAL]

    HISTORY_PROGRESS_TITLE = 'Progress'
    HISTORY_HIDDEN_TRACK_STORAGE = 'GSuite track storage'
    HISTORY_ERROR_TITLE = 'GSuite - files that failed manipulation'

    DEFAULT_SELECT = '--- Select ---'

    @staticmethod
    def getToolName():
        return "Expand primary tracks in GSuite with LD information"

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def getInputBoxNames():
        return [
            ('Select GSuite from history', 'gSuite'),
            ('Select reference genome', 'refGenome'),
            ('Select r<sup>2</sup> threshold', 'rsquare')
        ]

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxRefGenome(choices):
        return ['hg19', 'hg38']

    @staticmethod
    def getOptionsBoxRsquare(choices):
        return [LDGSuiteGeneratorTool.DEFAULT_SELECT, '1.0', '0.9', '0.8', '0.7', '0']

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
        from quick.webtools.clustering.CreateLDTrack import CreateLDTrack
        from quick.webtools.clustering.RsidMapper import RsidMapper

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        outGSuite = GSuite()
        errorGSuite = GSuite()

        progressViewer = ProgressViewer([('Manipulate tracks', gSuite.numTracks() + 24)],
                                        cls.extraGalaxyFn[cls.HISTORY_PROGRESS_TITLE])

        ldDict = CreateLDTrack.getLDDict(float(choices.rsquare))
        rsidDict = RsidMapper.createRsidMappingFromStaticFiles(progressViewer, choices.refGenome)
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
                    trackType='points',
                    attributes=track.attributes
                )

                trackFn = gSuiteTrack.path
                CreateLDTrack.parseFileIntoPointTrack(track.path, trackFn, ldDict, rsidDict)
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

        if choices.rsquare == LDGSuiteGeneratorTool.DEFAULT_SELECT:
            return 'Please select a threshold of r<sup>2<sup>'


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
        core.header("Convert GSuite to point track GSuite with LD information")
        core.paragraph("This tool takes in a GSuite, and generate a new GSuite of point tracks corresponding to each "
                       "of the original primary tracks. The point tracks of the new GSuite contain all variants that "
                       "are in LD with the original track elements.")
        core.paragraph("The r<sup>2</sup> "
                       "threshold sets the lower limit of LD correlation. r<sup>2</sup> >= 0.8 is the consensus "
                       "threshold for a significant LD correlation. Only variants with r<sup>2</sup> >= the chosen "
                       "threshold will be included in the point tracks.")
        core.divider()
        core.smallHeader("Primary track requirements")
        core.paragraph("The primary tracks need a column header of 'snps'. This column must contain rsids, which "
                       "represent the different track elements.")
        core.divider()
        core.smallHeader("Limitations")
        core.paragraph("The tool makes use of LD information from a static file within the HyperBrowser to expand the "
                       "primary tracks with LD information. This file was generated with all significant SNPs from the "
                       "GWAS Catalog (18660 from last extraction). If the GSuite tracks contain rsids that are not "
                       "present in the LD master file, the track element will remain in the track "
                       "as a single point without any "
                       "known LD variants")
        core.paragraph("The master LD file was created for the CEU population. Thus, this tool is not appropriate if "
                       "information of LD between SNPs of other populations is wanted.")
        core.divider()
        core.smallHeader("Example use")
        core.paragraph("Can be used to expand a GWAS Catalog GSuite with LD variants. One purpose of this expansion "
                       "could be to check overlap of LD variant loci against other genomic annotation tracks.")
        core.paragraph(
            "<b>NB:</b> Tool takes some time to run, as it reads in a full rsid-position mapping for all dbSNP rsids "
            "in the chosen reference genome.")
        return str(core)
