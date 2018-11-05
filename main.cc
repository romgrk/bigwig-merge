/*
 * main.c
 * Copyright (C) 2018 rgregoir <rgregoir@laurier>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdlib.h>
#include <string.h>
#include <map>
#include <cmath>

#include "./third_party/libBigWig/bigWig.h"


#define FOR_EACH_FILE(varname, block) for (size_t varname = 0; varname < filesLength; varname++) block


struct StringCompare
{
   bool operator()(char const *a, char const *b) const
   {
      return strcmp(a, b) < 0;
   }
};

void init() {
    // Initialize enough space to hold 128KiB (1<<17) of data at a time
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "Received an error in bwInit\n");
        exit(1);
    }
}

bigWigFile_t *openBigwigFile(char *filepath, const char *mode = "r") {
    auto file = bwOpen(filepath, NULL, mode);

    if (!file) {
        fprintf(stderr, "An error occured while opening %s\n", filepath);
        exit(2);
        return NULL;
    }

    return file;
}

int getChromLength(bigWigFile_t *file, const char * chrom) {
    auto chromList = file->cl;
    for (int i = 0; i < chromList->nKeys; i++) {
        if (strcmp(chromList->chrom[i], chrom) == 0)
            return chromList->len[i];
    }
    return -1;
}

double getChromMax(bigWigFile_t *file, char* chrom, uint32_t start, uint32_t end) {
    double* stats = bwStats(file, chrom, start, end, 1, bwStatsType::max);
    auto max = stats[0];
    free(stats);
    return max;
}

void printChromList(bigWigFile_t *file) {
    printf("Chroms:\n");

    // Print chrom list
    auto chromList = file->cl;
    for (int i = 0; i < chromList->nKeys; i++) {
        printf("\t%i:\t%s\t%i\n", i, chromList->chrom[i], chromList->len[i]);
    }

    printf("\n");
}

void printChromRange(bigWigFile_t *file, const char* chrom, uint32_t start, uint32_t end) {
    auto iter = bwOverlappingIntervalsIterator(file, (char *)chrom, start, end, 1);

    while (iter != NULL && iter->data) {

        printf("%i:%i-%i %p \n", iter->tid, iter->start, iter->end, iter->intervals);

        auto intervals = iter->intervals;
        if (intervals == NULL) {
            printf("\t(NULL) intervals \n" );
        }
        else {
            printf("\tLength: %i \n", intervals->l);
            for (uint32_t i = 0; i < intervals->l; i++) {
               printf("\t%i-%i: %f \n", intervals->start[i], intervals->end[i], intervals->value[i]);
            }
        }

        auto entries = iter->entries;
        if (entries == NULL) {
            printf("\t(NULL) entries \n" );
        }
        else {
            printf("\tLength: %i \n", entries->l);
            // for (uint32_t i = 0; i < entries->l; i++) {
                // printf("\t%i-%i: %f \n", entries->start[i], entries->end[i], entries->value[i]);
            // }
        }

        iter = bwIteratorNext(iter);
    }

    bwIteratorDestroy(iter);
}

struct ChromDetails {
    bool *hasFile;
    uint32_t length;
};
typedef struct ChromDetails ChromDetails;
typedef std::map<char*, ChromDetails, StringCompare> ChromDetailsMap;

ChromDetailsMap getChromDetails(bigWigFile_t **files, size_t length) {
    ChromDetailsMap chromDetails;

    auto file1 = files[0];

    auto chromList = file1->cl;
    for (int i = 0; i < chromList->nKeys; i++) {

        auto name = chromList->chrom[i];
        auto hasFile = new bool[length]();
        hasFile[0] = true;
        auto details = (ChromDetails){
            .hasFile = hasFile,
            .length = chromList->len[i]
        };
        chromDetails[name] = details;
    }

    for (size_t n = 0; n < length; n++) {
        auto file = files[n];

        auto chromList = file->cl;
        for (int i = 0; i < chromList->nKeys; i++) {

            auto name = chromList->chrom[i];

            if (chromDetails.find(name) == chromDetails.end()) {
                auto name = chromList->chrom[i];
                auto hasFile = new bool[length]();
                hasFile[n] = true;
                auto details = (ChromDetails){
                    .hasFile = hasFile,
                    .length = chromList->len[i]
                };
                chromDetails[name] = details;
            } else {
                auto details = chromDetails[name];
                details.hasFile[n] = true;
                chromDetails[name] = details;
            }
        }
    }


    return chromDetails;
}

void createChromList(bigWigFile_t *file, ChromDetailsMap chromDetails) {
    size_t size = chromDetails.size();

    char **chroms = new char*[size];
    uint32_t *lengths = new uint32_t[size];

    auto iterator = chromDetails.begin();
    int index = 0;
    while (iterator != chromDetails.end()) {
        auto name = iterator->first;
        auto details = iterator->second;

        chroms[index] = name;
        lengths[index] = details.length;

        iterator++;
        index++;
    }

    file->cl = bwCreateChromList(chroms, lengths, (int64_t)size);
}

int main(int argc, char *argv[]) {

    if (argc < 4) {
        fprintf(stderr, "Usage: %s file1.bw file2.bw .. fileN.bw output.bw\n", argv[0]);
        return 1;
    }

    init();

    // Open the local/remote file

    size_t filesLength = argc - 2;
    bigWigFile_t **files = new bigWigFile_t*[filesLength];

    FOR_EACH_FILE(n, { files[n] = openBigwigFile(argv[1]); })

    auto output = openBigwigFile(argv[3], "w");

    /* printChromList(file1);
     * printChromList(file2);
     * printf("\n"); */

    ChromDetailsMap chromDetails = getChromDetails(files, filesLength);
    ChromDetailsMap::iterator iterator;

    double *factors = NULL;
    double totalMaximum = 0;


    if (bwCreateHdr(output, 10)) goto error;

    createChromList(output, chromDetails);
    if (!output->cl) goto error;

    // Write the header
    if (bwWriteHdr(output)) goto error;


    /*
     * Calculate factor for values of each file
     */

    factors = new double[filesLength];
    FOR_EACH_FILE(n, {
        factors[n] = files[n]->hdr->maxVal;
    })
    printf("\nMaximums: ");
    for (size_t i = 0; i < filesLength; i++) {
        totalMaximum += factors[i];
        printf("%f, ", factors[i]);
    }
    printf("\nTotalMaximum: %f", totalMaximum);
    printf("\nFactors: ");
    for (size_t i = 0; i < filesLength; i++) {
        factors[i] = factors[i] / totalMaximum;
        printf("%f, ", factors[i]);
    }
    printf("\n");


    iterator = chromDetails.begin();
    while(iterator != chromDetails.end()) {
        auto chrom = iterator->first;
        auto details = iterator->second;

        printf("\n%s\n", chrom);
        printf("\t%s: file1=%s file2=%s length=%i\n",
                chrom,
                details.hasFile[0] ? "true" : "false",
                details.hasFile[1] ? "true" : "false",
                details.length);

        bwOverlappingIntervals_t **intervals = new bwOverlappingIntervals_t*[filesLength]();

        FOR_EACH_FILE(n, {
            if (details.hasFile[n]) {
                intervals[n] = bwGetValues(files[n], chrom, 0, details.length, 1);
            }
        })

        uint32_t *starts = new uint32_t[details.length];
        uint32_t *ends   = new uint32_t[details.length];
        float *values    = new float[details.length];

        // Number of values
        uint32_t count = 0;

        for (uint32_t i = 0; i < details.length; i++) {
            bool isAllNaN = true;
            float value = 0;

            FOR_EACH_FILE(n, {
                if (std::isnan(intervals[n]->value[i]))
                    continue;

                isAllNaN = false;

                value += intervals[n]->value[i] * factors[n];
            })

            if (isAllNaN)
                continue;

            starts[count] = i;
            ends[count] = i + 1;
            values[count] = value;

            count++;
        }

        printf("\tAdding %i intervals to %s\n", count, chrom);

        if (bwAddIntervalsToChrom(output, chrom, starts, ends, values, count)) goto error_add;

        delete[] intervals;
        delete[] starts;
        delete[] ends;
        delete[] values;
        iterator++;

        continue;

    error_add:
        fprintf(stderr, "\x1b[90mbwAddIntervalToChrom\x1b[0m\n");
        goto error;
    }


    // chr6: 52860296, 52860836
    // chr6: 52848630, 53009526
    // printChromRange(files[0], "chr6", 52848630, 53009526);


    // Get values in a range (0-based, half open) without NAs
    /* intervals = bwGetValues(files[0], "chr1", 60000, 65000, 0);
     * bwDestroyOverlappingIntervals(intervals); // Free allocated memory */

    //Get values in a range (0-based, half open) with NAs
    /* intervals = bwGetValues(files[0], "chr1", 10000000, 10000100, 1);
     * bwDestroyOverlappingIntervals(intervals); // Free allocated memory */

    //Get the full intervals that overlap
    /* intervals = bwGetOverlappingIntervals(files[0], "chr1", 10000000, 10000100);
     * bwDestroyOverlappingIntervals(intervals); */

    //Get an example statistic - standard deviation
    //We want ~4 bins in the range
    /* double *stats = bwStats(files[0], "chr1", 10000000, 10000100, 4, dev);
     * if(stats) {
     *     printf("chr1:10000000-10000100 std. dev.: %f %f %f %f\n", stats[0], stats[1], stats[2], stats[3]);
     *     free(stats);
     * } */


    FOR_EACH_FILE(n, { bwClose(files[n]); })
    bwClose(output);
    bwCleanup();
    delete[] files;
    delete[] factors;

    return 0;

error:
    FOR_EACH_FILE(n, { bwClose(files[n]); })
    bwClose(output);
    bwCleanup();
    delete[] files;
    delete[] factors;

    fprintf(stderr, "\x1b[90mReceived error\x1b[0m\n");

    return 1;
}
