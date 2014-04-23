#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append(sys.path[0] + "/res")
import re
import urllib.request
from Bio.Blast import NCBIWWW, NCBIXML
from gi.repository import Gtk, GObject
import threading


class Iface(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="BLASTats")
        self.set_resizable(False)
        self.set_default_size(800, 600)
        self.set_position(Gtk.WindowPosition.CENTER)
        self.connect("key-press-event", self.on_key_press)

        grid = Gtk.Grid()
        grid.set_row_spacing(10)
        self.add(grid)

        grid_sequence = Gtk.Grid()
        label_sequence = Gtk.Label('Sequence:')
        label_sequence.set_alignment(0, 0.5)
        self.entry_sequence = Gtk.Entry()
        self.entry_sequence.set_width_chars(46)
        grid_sequence.add(label_sequence)
        grid_sequence.attach(self.entry_sequence, 1, 0, 1, 1)

        frame_options = Gtk.Frame(label="Options")
        frame_options.set_hexpand(True)
        grid_options = Gtk.Grid()
        grid_options.set_hexpand(True)
        label_idthresh = Gtk.Label('Identity threshold: ')
        label_idthresh.set_alignment(0, 0.5)
        label_covthresh = Gtk.Label('Coverage threshold: ')
        label_pc1 = Gtk.Label("%                                  ")
        label_pc2 = Gtk.Label("%                                  ")
        self.entry_idthresh = Gtk.Entry()
        self.entry_idthresh.set_max_length(3)
        self.entry_idthresh.set_width_chars(3)
        self.entry_idthresh.set_placeholder_text("80")
        self.entry_covthresh = Gtk.Entry()
        self.entry_covthresh.set_max_length(3)
        self.entry_covthresh.set_width_chars(3)
        self.entry_covthresh.set_placeholder_text("95")
        self.check_verbose = Gtk.CheckButton("Verbose")
        self.check_verbose.set_active(True)
        self.check_list = Gtk.CheckButton("List all species")
        grid_options.add(label_idthresh)
        grid_options.attach(label_covthresh, 0, 1, 1, 1)
        grid_options.attach(self.entry_idthresh, 1, 0, 1, 1)
        grid_options.attach(self.entry_covthresh, 1, 1, 1, 1)
        grid_options.attach(label_pc1, 2, 0, 1, 1)
        grid_options.attach(label_pc2, 2, 1, 1, 1)
        grid_options.attach(self.check_verbose, 3, 0, 1, 1)
        grid_options.attach(self.check_list, 3, 1, 1, 1)
        frame_options.add(grid_options)

        grid_buttons = Gtk.Grid()
        button_clear = Gtk.Button("Clear output")
        button_analyse = Gtk.Button("Analyse")
        button_analyse.set_size_request(150, -1)
        button_blast = Gtk.Button("Blast & Analyse")
        button_blast.connect("clicked", self.on_blast_click)
        button_blast.set_size_request(150, -1)
        button_analyse.connect("clicked", self.on_analyse_click)
        button_clear.connect("clicked", self.clear_output)
        button_clear.set_size_request(150, -1)
        grid_buttons.attach(button_clear, 0, 0, 1, 1)
        grid_buttons.attach(button_analyse, 1, 0, 1, 1)
        grid_buttons.attach(button_blast, 2, 0, 1, 1)

        frame_output = Gtk.Frame(label="Output")
        self.txtview_output = Gtk.TextView()
        self.txtview_output.set_editable(False)
        scroll_output = Gtk.ScrolledWindow()
        scroll_output.set_min_content_height(400)
        scroll_output.set_min_content_width(400)
        scroll_output.add(self.txtview_output)
        frame_output.add(scroll_output)

        grid_footer = Gtk.Grid()
        button_help = Gtk.Button(self, stock="gtk-help")
        button_help.connect("clicked", self.help_)
        button_quit = Gtk.Button(self, stock="gtk-quit")
        button_quit.connect("clicked", Gtk.main_quit)
        grid_footer.add(button_help)
        grid_footer.attach(button_quit, 1, 0, 1, 1)

        grid.add(grid_sequence)
        grid.attach(frame_options, 0, 1, 1, 1)
        grid.attach(grid_buttons, 0, 2, 1, 1)
        grid.attach(frame_output, 0, 3, 1, 1)
        grid.attach(grid_footer, 0, 4, 1, 1)

    def print_(self, txt):
        buf = self.txtview_output.get_buffer()
        buf.insert(buf.get_end_iter(), txt + '\n')

    def help_(self, *args):
        self.clear_output()
        self.print_("Description\n----------------------\n"
                    "BLASTats helps you see how your protein of interest is distributed among B. cereus group.\n\n"
                    "Parameters\n---------------------\n"
                    "- Identity threshold:\nYou can set manually the identity threshold required between your query and"
                    " a BLAST hit.\nOnly hits above the threshold will be considered as potential functional homologues"
                    " of your query protein.\n"
                    "- Coverage threshold:\nYou can also set manually the query coverage threshold required for a"
                    " BLAST hit to be considered.\n\n"
                    "References\n----------------------\n"
                    "For more information about functional homology infering, see:\n"
                    "Rost, B., Liu, J., Nair, R., Wrzeszczynski, K. O., & Ofran, Y. (2003). Automatic prediction of"
                    " protein function. Cellular and Molecular Life Sciences : CMLS, 60(12), 2637–50.")

    def clear_output(self, *args):
        buf = self.txtview_output.get_buffer()
        buf.set_text("")

    def fetch_arguments(self, check_seq=False):
        kwargs = {}

        if check_seq:
            seq = self.entry_sequence.get_text()
            anti_aminoacids = re.compile('[^AaBbCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtUuVvWwYyZzXx\*\-]')
            if len(anti_aminoacids.findall(seq)) == 0:
                kwargs['seq'] = seq
            else:
                self.print_("Error:\nPlease provide a valid protein sequence\n"
                            "Allowed characters: ABCDEFGHIKLMNPQRSTUVWXYZ*-\n")

        if self.check_verbose.get_active():
            kwargs['verbose'] = True

        if self.check_list.get_active():
            kwargs['details'] = True

        if self.entry_idthresh.get_text() != "":
            try:
                if 100 >= int(self.entry_idthresh.get_text()) >= 0:
                    kwargs['identity_threshold'] = int(self.entry_idthresh.get_text())/100
                else:
                    self.print_("Error:\nIdentity threshold must be between 0 and 100")
            except ValueError:
                self.print_("Error:\nPlease provide a valid identity threshold")

        if self.entry_covthresh.get_text() != "":
            try:
                if 100 >= int(self.entry_covthresh.get_text()) >= 0:
                    kwargs['query_cover_threshold'] = int(self.entry_covthresh.get_text())/100
                else:
                    self.print_("Error:\nQuery coverage threshold must be between 0 and 100")
            except ValueError:
                self.parent.print_("Error:\nPlease provide a query coverage threshold")

        return kwargs

    def on_blast_click(self, *args):
        kwargs = self.fetch_arguments(check_seq=True)
        if kwargs['seq'] != "":
            compute = Compute(self, **kwargs)
            if self.check_verbose.get_active():
                self.print_("BLASTing sequence at NCBI...")
            threading.Thread(target=compute.blast).start()
        else:
            self.print_("Error\nPlease enter a protein sequence to blast\n")

    def on_analyse_click(self, *args):
        kwargs = self.fetch_arguments()
        compute = Compute(self, **kwargs)
        threading.Thread(target=compute.analyse).start()

    @staticmethod
    def on_key_press(widget, event):
        if event.keyval == 65293 or event.keyval == 65421:
            pass
        if event.keyval == 65307:
            Gtk.main_quit()


class Compute(object):
    def __init__(self, parent, seq="", verbose=False, details=False, identity_threshold=0.8,
                 query_cover_threshold=0.95):
        self.parent = parent
        self.seq = seq
        self.verbose = verbose
        self.details = details
        self.identity_threshold = identity_threshold
        self.query_cover_threshold = query_cover_threshold

    @staticmethod
    def fetch_organism(string):
        """
        Returns a list containing organisms names contained in BLAST match title.
        Names are filtered to ignore names < 2 words (like "Bacillus cereus"), names containing "group" and names
        repeated several times in a single BLAST match title.
        """

        pattern = re.compile('\[(.*?)\]')
        matchlist = []
        for match in pattern.findall(string):
            if match not in matchlist and "group" not in match and len(match.split()) > 2:
                matchlist.append(match)

        return matchlist

    @staticmethod
    def fetch_genomes_quantity():
        """
        Returns a dictionary containing the number of sequenced genomes available at NCBI for Bc, Bt, Ba and Bs.
        """

        try:

            genomes_quantity = {}
            for organism in ['bacillus+cereus', 'bacillus+thuringiensis', 'bacillus+anthracis', 'bacillus+subtilis']:
                ans_bc = urllib.request.urlopen("http://www.ncbi.nlm.nih.gov/genome/?term={}".format(organism))
                page_bc = ans_bc.read()
                # Because ans.read() returns a byte string while regex findall takes unicode strings
                page_bc = page_bc.decode()

                # What is retrieved from NCBI's page here are 3 occurences of numbers between brackets being
                # respectively the number of genomes, scaffolds/contigs and SRA/traces
                pattern = re.compile('\[(\d+)\]')
                grep = pattern.findall(page_bc)
                # Convert strings list to ints list
                for i in range(len(grep)):
                    grep[i] = int(grep[i])

                if organism == 'bacillus+cereus':
                    genomes_quantity['bc'] = sum(grep)
                elif organism == 'bacillus+thuringiensis':
                    genomes_quantity['bt'] = sum(grep)
                elif organism == 'bacillus+anthracis':
                    genomes_quantity['ba'] = sum(grep)
                elif organism == 'bacillus+subtilis':
                    genomes_quantity['bs'] = sum(grep)

            return genomes_quantity

        except urllib.request.URLError:
            print("\nError:\nCould not reach NCBI's website.\nPlease check your internet connection and retry.\n")
            sys.exit(0)

    def blast(self):
        """
        BLASTs protein sequence at NCBI and writes the result in a "blast_results.xml" file in the execution directory.
        """

        try:
            query = NCBIWWW.qblast('blastp', 'nr', self.seq, hitlist_size=500)
        except urllib.error.URLError:
            self.parent.print_("\nError:\nCould not reach NCBI's website.\n"
                               "Please check your internet connection and retry.\n")
        else:
            of_ = open('blast_results.xml', 'w')
            of_.write(query.read())
            of_.close()
            query.close()
            self.analyse()

    def analyse(self):

        if self.verbose:
            GObject.idle_add(self.parent.print_, "Setting identity threshold to {}".format(self.identity_threshold))
            GObject.idle_add(self.parent.print_, "Setting query coverage threshold to {}".
                             format(self.query_cover_threshold))

        # Lists containing hits for each member
        total, bacilli, bc, bt, ba, bs, sp = [], [], [], [], [], [], []
        # Dictionary containing number of sequenced genomes for each member
        if self.verbose:
            GObject.idle_add(self.parent.print_, "Retrieving genomes quantities at NCBI...")
        genomes = self.fetch_genomes_quantity()

        if self.verbose:
            GObject.idle_add(self.parent.print_, "Parsing BLAST results...")

        try:
            results_handle = open('blast_results.xml')
            results = NCBIXML.read(results_handle)
        except FileNotFoundError:
            self.parent.print_("\nError:\nResults file could not be found. Please BLAST a sequence and retry.\n")
        except ValueError:
            self.parent.print_("\nError:\nCould not read your XML results file (is it empty?). "
                               "Please check it and retry\n")
        else:

            for result in results.alignments:
                hsp = result.hsps[0]
                query_cover = (len(hsp.sbjct) - hsp.sbjct_start) / results.query_letters
                identity = hsp.identities / len(hsp.sbjct)
                if query_cover > self.query_cover_threshold and identity > self.identity_threshold:
                    organism_list = self.fetch_organism(result.title)
                    for organism in organism_list:
                        if organism not in total:
                            total.append(organism)
                            if "Bacillus" in organism:
                                bacilli.append(organism)
                                if "cereus" in organism:
                                    bc.append(organism)
                                elif "thuringiensis" in organism:
                                    bt.append(organism)
                                elif "anthracis" in organism:
                                    ba.append(organism)
                                elif "subtilis" in organism:
                                    bs.append(organism)
                                else:
                                    sp.append(organism)

            GObject.idle_add(self.parent.print_, '==============================')
            GObject.idle_add(self.parent.print_, "Total: {}".format(len(total)))
            GObject.idle_add(self.parent.print_, "|-Bacilli: {}".format(len(bacilli)))
            GObject.idle_add(self.parent.print_, "|--- cereus: {0} ({1}%)".
                             format(len(bc), int(len(bc) / genomes['bc'] * 100)))
            GObject.idle_add(self.parent.print_, "|--- thuringiensis: {0} ({1}%)".
                             format(len(bt), int(len(bt) / genomes['bt'] * 100)))
            GObject.idle_add(self.parent.print_, "|--- anthracis: {0} ({1}%)".
                             format(len(ba), int(len(ba) / genomes['ba'] * 100)))
            GObject.idle_add(self.parent.print_, "|--- subtilis: {0} ({1}%)".
                             format(len(bs), int(len(bs) / genomes['bs'] * 100)))
            GObject.idle_add(self.parent.print_, "|--- Other Bacilli: {}".format(len(sp)))
            GObject.idle_add(self.parent.print_, "|-Other non-Bacilli: {}".format(len(total)-len(bacilli)))
            GObject.idle_add(self.parent.print_, '==============================')

            if self.details:
                bc.sort()
                bt.sort()
                ba.sort()
                bs.sort()
                sp.sort()
                GObject.idle_add(self.parent.print_, "\nBacillus cereus:\n-----------------------")
                for bac in bc:
                    GObject.idle_add(self.parent.print_, bac)
                GObject.idle_add(self.parent.print_, "\nBacillus thuringiensis:\n-----------------------")
                for bac in bt:
                    GObject.idle_add(self.parent.print_, bac)
                GObject.idle_add(self.parent.print_, "\nBacillus anthracis:\n-----------------------")
                for bac in ba:
                    GObject.idle_add(self.parent.print_, bac)
                GObject.idle_add(self.parent.print_, "\nBacillus subtilis:\n-----------------------")
                for bac in bs:
                    GObject.idle_add(self.parent.print_, bac)
                GObject.idle_add(self.parent.print_, "\nOther Bacilli:\n-----------------------")
                for bac in sp:
                    GObject.idle_add(self.parent.print_, bac)


def main():
    iface = Iface()
    iface.connect("delete-event", Gtk.main_quit)
    iface.show_all()
    Gtk.main()


if __name__ == '__main__':
    main()
