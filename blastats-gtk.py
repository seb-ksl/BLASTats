#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  pystrains.py
#
#  Copyright 2013-2014 Sébastien Gélis-Jeanvoine <sebastien@gelis.ch>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.
#

# FixMe: cannot fetch when only 1 genome available (cytotoxicus & toyonensis)
#
# ToDo: Pickle ethresh
# ToDo: add chromosomes to analysis (not only complete genomes), but how to build entrez query on NCBI's BLAST ?
# ToDo: add possibility to filter by e-value.
# ToDo: refine name fetching. For example, "Clostridium difficile complete genome, strain M120"
#       will be treated as "Clostridium difficile"
# ToDo: decrease number of hsps per hit in BLAST arguments
# ToDo: enhance GUI so user can choose if aligning HSPs or retrieved CDSs
# ToDo: capture memory overflow when parsing too big BLAST results
# ToDo: conditional imports with try

import sys
sys.path.append("{}/res".format(sys.path[0]))

# from Bio import Phylo
from Bio import Entrez, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.Applications._Fasttree import FastTreeCommandline

from gi.repository import Gtk, GObject, Pango
import pickle
import queue
import re
import subprocess
import threading
import urllib.request

# GLOBALS
GENOMES_URL = "http://www.ncbi.nlm.nih.gov/genome/genomes/{}"
GENOMES_REGEX = "Complete \[(\d+)\]"
# Settings
Entrez.email = "sgelis@jouy.inra.fr"


class Iface(Gtk.Window):
    """
    GTK interface class.
    """
    def __init__(self):
        Gtk.Window.__init__(self, title="BLASTats")
        self.set_resizable(False)
        self.set_default_size(800, 600)
        self.set_position(Gtk.WindowPosition.CENTER)
        self.connect("key-press-event", self.on_key_press)
        self.connect("delete-event", self.quit)

        grid = Gtk.Grid()
        grid.set_row_spacing(10)
        self.add(grid)

        grid_sequence = Gtk.Grid()
        label_sequence = Gtk.Label("Sequence:")
        label_sequence.set_alignment(0, 0.5)
        self.entry_sequence = Gtk.Entry()
        self.entry_sequence.set_width_chars(46)
        grid_sequence.add(label_sequence)
        grid_sequence.attach(self.entry_sequence, 1, 0, 1, 1)

        grid_organisms = Gtk.Grid()
        label_organism = Gtk.Label("Organism: ")
        self.entry_organism = Gtk.Entry()
        self.entry_organism.set_placeholder_text("Organism name")
        self.entry_organism.connect("key-press-event", self.on_organism_add)
        scroll_organism = Gtk.ScrolledWindow()
        scroll_organism.set_min_content_height(60)
        scroll_organism.set_min_content_width(200)
        self.liststore_organism = Gtk.ListStore(str)
        renderer_organism = Gtk.CellRendererText()
        self.treeview_organism = Gtk.TreeView(model=self.liststore_organism)
        self.treeview_organism.append_column(Gtk.TreeViewColumn("Organism", renderer_organism, text=0))
        self.treeview_organism.set_headers_visible(False)
        self.treeview_organism.connect("button-press-event", self.on_treeview_click)
        scroll_organism.add(self.treeview_organism)
        grid_organisms.add(label_organism)
        grid_organisms.attach(self.entry_organism, 1, 0, 1, 1)
        grid_organisms.attach(scroll_organism, 2, 0, 1, 1)

        frame_options = Gtk.Frame(label="Options")
        grid_options = Gtk.Grid()
        label_idthresh = Gtk.Label("Identity threshold: ")
        label_idthresh.set_alignment(0, 0.5)
        label_covthresh = Gtk.Label("Coverage threshold: ")
        label_ethresh = Gtk.Label("E-value threshold: ")
        label_pc1 = Gtk.Label("%         ")
        label_pc2 = Gtk.Label("%         ")
        self.entry_idthresh = Gtk.Entry()
        self.entry_idthresh.set_max_length(3)
        self.entry_idthresh.set_width_chars(3)
        self.entry_idthresh.set_max_width_chars(3)
        self.entry_idthresh.set_halign(Gtk.Align.START)
        self.entry_idthresh.set_size_request(20, -1)
        self.entry_idthresh.set_placeholder_text("80")
        self.entry_covthresh = Gtk.Entry()
        self.entry_covthresh.set_max_length(3)
        self.entry_covthresh.set_width_chars(3)
        self.entry_covthresh.set_max_width_chars(3)
        self.entry_covthresh.set_halign(Gtk.Align.START)
        self.entry_covthresh.set_placeholder_text("95")
        self.entry_ethresh = Gtk.Entry()
        self.entry_ethresh.set_width_chars(5)
        self.entry_ethresh.set_max_width_chars(5)
        self.entry_ethresh.set_halign(Gtk.Align.START)
        self.entry_ethresh.set_placeholder_text("1e-5")
        self.check_verbose = Gtk.CheckButton("Verbose")
        self.check_list = Gtk.CheckButton("List all species")
        self.check_fasta = Gtk.CheckButton("Save results in FASTA")
        self.check_align_tree = Gtk.CheckButton("Align & Tree")
        self.check_align_tree.connect("toggled", self.on_align_tree_click)
        grid_options.add(label_idthresh)
        grid_options.attach(label_covthresh, 0, 1, 1, 1)
        grid_options.attach(label_ethresh, 0, 2, 1, 1)
        grid_options.attach(self.entry_idthresh, 1, 0, 1, 1)
        grid_options.attach(self.entry_covthresh, 1, 1, 1, 1)
        grid_options.attach(self.entry_ethresh, 1, 2, 1, 1)
        grid_options.attach(label_pc1, 2, 0, 1, 1)
        grid_options.attach(label_pc2, 2, 1, 1, 1)
        grid_options.attach(self.check_verbose, 3, 0, 1, 1)
        grid_options.attach(self.check_list, 3, 1, 1, 1)
        grid_options.attach(self.check_fasta, 3, 2, 1, 1)
        grid_options.attach(self.check_align_tree, 3, 3, 1, 1)
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
        self.txtview_output.modify_font(Pango.FontDescription("mono"))
        scroll_output = Gtk.ScrolledWindow()
        scroll_output.set_min_content_height(400)
        scroll_output.set_min_content_width(400)
        scroll_output.add(self.txtview_output)
        frame_output.add(scroll_output)

        grid_footer = Gtk.Grid()
        button_help = Gtk.Button(self, stock="gtk-help")
        button_help.connect("clicked", self.help_)
        button_quit = Gtk.Button(self, stock="gtk-quit")
        button_quit.connect("clicked", self.quit)
        grid_footer.add(button_help)
        grid_footer.attach(button_quit, 1, 0, 1, 1)

        grid.add(grid_sequence)
        grid.attach(grid_organisms, 0, 1, 1, 1)
        grid.attach(frame_options, 0, 2, 1, 1)
        grid.attach(grid_buttons, 0, 3, 1, 1)
        grid.attach(frame_output, 0, 4, 1, 1)
        grid.attach(grid_footer, 0, 5, 1, 1)

        # Although list_organisms can seem redundant with liststore_organism, it allows checking for duplicates
        self.list_organisms = []
        self.read_settings()

    def read_settings(self):
        try:
            fhandle = open("settings", "rb")
        except FileNotFoundError:
            print("Settings not found")
        else:
            settings = pickle.load(fhandle)
            self.check_verbose.set_active(settings["verbose"])
            self.check_list.set_active(settings["list"])
            self.check_fasta.set_active(settings["fasta"])
            self.check_align_tree.set_active(settings["align_tree"])
            self.entry_idthresh.set_text(settings["idthresh"])
            self.entry_covthresh.set_text(settings["covthresh"])
            self.entry_ethresh.set_text(settings["ethresh"])
            for organism in settings["organisms"]:
                self.list_organisms.append(organism)
                self.liststore_organism.append((organism,))

    def print_(self, txt):
        buf = self.txtview_output.get_buffer()
        buf.insert(buf.get_end_iter(), "{}\n".format(str(txt)))

    def help_(self, *args):
        self.clear_output()
        self.print_("For help, see http://www.gelis.ch/programs/blastats/")

    def clear_output(self, *args):
        buf = self.txtview_output.get_buffer()
        buf.set_text("")

    def fetch_arguments(self, check_seq=False):
        kwargs = {}

        if check_seq:
            seq = self.entry_sequence.get_text().replace("\n", "")
            anti_aminoacids = re.compile("[^AaBbCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtUuVvWwYyZzXx\*\-]")
            if len(anti_aminoacids.findall(seq)) == 0:
                kwargs["seq"] = seq
            else:
                self.print_("Error:\nPlease provide a valid protein sequence\n"
                            "Allowed characters: ABCDEFGHIKLMNPQRSTUVWXYZ*-\n")

        organisms = []
        for row in self.liststore_organism:
            organisms.append(row[0].replace(" ", "+"))
        organisms.sort()
        kwargs["organisms"] = organisms

        kwargs["verbose"] = self.check_verbose.get_active()
        kwargs["details"] = self.check_list.get_active()
        kwargs["fasta"] = self.check_fasta.get_active()
        kwargs["align_tree"] = self.check_align_tree.get_active()

        if self.entry_idthresh.get_text() != "":
            try:
                if 100 >= int(self.entry_idthresh.get_text()) >= 0:
                    kwargs["identity_threshold"] = int(self.entry_idthresh.get_text())/100
                else:
                    self.print_("Error:\nIdentity threshold must be between 0 and 100")
            except ValueError:
                self.print_("Error:\nPlease provide a valid identity threshold")

        if self.entry_covthresh.get_text() != "":
            try:
                if 100 >= int(self.entry_covthresh.get_text()) >= 0:
                    kwargs["query_cover_threshold"] = int(self.entry_covthresh.get_text()) / 100
                else:
                    self.print_("Error:\nQuery coverage threshold must be between 0 and 100")
            except ValueError:
                self.parent.print_("Error:\nPlease provide a query coverage threshold")

        if self.entry_ethresh.get_text() != "":
            try:
                if 1 >= float(self.entry_ethresh.get_text()) >= 0:
                    kwargs["e_threshold"] = float(self.entry_ethresh.get_text())
                else:
                    self.print_("Error:\nE-value threshold must be between 0 and 1")
            except ValueError:
                self.parent.print_("Error:\nPlease provide an E-value threshold")

        return kwargs

    def on_organism_add(self, widget, event):
        if event.keyval == 65293 or event.keyval == 65421:
            organism = self.entry_organism.get_text().capitalize()
            if len(organism.split()) == 2:
                if organism not in self.list_organisms:
                    self.liststore_organism.append((organism,))
                else:
                    self.print_("Error: organism already in list.")
            else:
                self.print_("Error: organism name must be 2 words long.")

    def on_treeview_click(self, widget, event):
        if event.button == 3:
            # Determine what row is under the cursor. path[0] is the row number, path[1] is the column,
            # path[3] is cell(x) and path[4] is cell(y)
            path = self.treeview_organism.get_path_at_pos(event.x, event.y)
            if path is not None:
                # Set keyboard focus to treeview
                self.treeview_organism.grab_focus()
                # Set keyboard focus to the right row and column
                self.treeview_organism.set_cursor(path[0], path[1], 0)
                model, row = self.treeview_organism.get_selection().get_selected()
                self.liststore_organism.remove(row)

    def on_align_tree_click(self, *args):
        if self.check_align_tree.get_active():
            self.check_fasta.set_active(True)
            self.check_fasta.set_sensitive(False)
        else:
            self.check_fasta.set_sensitive(True)

    def on_blast_click(self, *args):
        kwargs = self.fetch_arguments(check_seq=True)
        if kwargs.get("seq", "") != "":
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

    def on_key_press(self, widget, event):
        if event.keyval == 65293 or event.keyval == 65421:
            pass
        if event.keyval == 65307:
            self.quit()

    def quit(self, *args):
        try:
            fhandle = open("settings", "wb")
        except PermissionError:
            print("Not allowed to write settings.")
        else:
            settings = {}
            settings["verbose"] = self.check_verbose.get_active()
            settings["list"] = self.check_list.get_active()
            settings["fasta"] = self.check_fasta.get_active()
            settings["align_tree"] = self.check_align_tree.get_active()
            settings["idthresh"] = self.entry_idthresh.get_text()
            settings["covthresh"] = self.entry_covthresh.get_text()
            settings["ethresh"] = self.entry_ethresh.get_text()
            settings["organisms"] = []
            for organism in self.liststore_organism:
                settings["organisms"].append(organism[0])
            pickle.dump(settings, fhandle)
        finally:
            fhandle.close()
            Gtk.main_quit()


class Compute(object):
    """
    Computing object.
    The compute object takes parameters given by the interface and:
    1) -Optional- BLASTs the protein sequence
    2) Fetches the number of sequenced genomes at NCBI's database
    3) Parses the BLAST's XML output to
        3.1) find the organisms in each sequence match
        3.2) build a tree of organisms given their genus and specie
    4) Outputs a tree of species of interest (inputted by interface), computing for each specie the abundance of the
       query sequence (given by the number of species listed in the matches tree under
       [genus_of_interest][specie_of_interest]) and the fraction of sequenced organisms it represents.
    """
    def __init__(self, parent, organisms, seq="", verbose=False, details=False, fasta=False, align_tree=False,
                 identity_threshold=0.8, query_cover_threshold=0.95, e_threshold=1e-5):
        self.parent = parent
        self.organisms_of_interest = organisms
        self.seq = seq
        self.verbose = verbose
        self.details = details
        self.fasta = fasta
        self.align_tree = align_tree
        self.identity_threshold = identity_threshold
        self.query_cover_threshold = query_cover_threshold
        self.e_threshold = e_threshold

    def get_url(self, organism, q):
        """
        Threaded method that gets a URL and puts its content in a threaded queue.
        """

        entrez_handle = Entrez.esearch(db="genome", term="{}[Organism]".format(organism))
        try:
            organism_id = Entrez.read(entrez_handle)["IdList"][0]
        except IndexError:
            q.put((organism, None))
        else:
            url = GENOMES_URL.format(organism_id)

        try:
            # Decode because .read() returns a byte string while re.findall() takes unicode strings
            page = urllib.request.urlopen(url).read().decode()
        except urllib.request.URLError:
            GObject.idle_add(self.parent.print_, "Error:\nCould not reach NCBI's website.\n"
                                                 "Please check your internet connection and retry.")
            q.put((organism, None))
        else:
            q.put((organism, page))

    def fetch_genomes_quantity(self):
        """
        Returns a dictionary containing the number of sequenced genomes available at NCBI for species of interest.
        URL requests are threaded via get_url(), and the results queue is read after all threads are done.
        """

        #genomes = {}
        genomes = {"Bacillus+toyonensis": 1, "Bacillus+bombysepticus": 1, "Bacillus+cytotoxicus": 1}
        q = queue.Queue()
        threads = []

        for organism in self.organisms_of_interest:
            threads.append(threading.Thread(target=self.get_url, args=(organism, q)))

        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        regex_genomes = re.compile(GENOMES_REGEX)

        while not q.empty():
            organism, page_organism = q.get()
            if isinstance(page_organism, type(None)) or "Organism Overview" not in page_organism:
                GObject.idle_add(self.parent.print_, "Could not find '{}' in NCBI's database. "
                                                     "Dropping it from analysis".format(organism.replace("+", " ")))
                self.organisms_of_interest.remove(organism)
            else:
                try:
                    genomes[organism] = int(regex_genomes.findall(page_organism)[0])
                except IndexError:
                    GObject.idle_add(self.parent.print_, "Error fetching genomes quantity for {}. Dropping it from "
                                                         "analysis".format(organism.replace("+", " ")))
                    if organism != "Bacillus+toyonensis" and organism != "Bacillus+bombysepticus" and organism != "Bacillus+cytotoxicus":
                        self.organisms_of_interest.remove(organism)
        return genomes

    @staticmethod
    def fetch_organism(string):
        """
        Returns a list containing organisms names contained in BLAST match title.
        Names are filtered to ignore names < 2 words (like "Bacillus cereus"), names containing "group" and names
        repeated several times in a single BLAST match title.
        """

        # pattern = re.compile("\[(.*?)\]")
        # matchlist = set()
        # for match in pattern.findall(string):
        #     if "group" not in match and len(match.split()) > 2:
        #         matchlist.add(match)

        organism = string.split("|")[-1].split(",")[0].lstrip()
        if "genom" in organism:
            organism = organism[:organism.find("genom")].rstrip()
        if "complete" in organism:
            organism = organism[:organism.find("complete")].rstrip()
        if "chromosome" in organism:
            organism = organism[:organism.find("chromosome")].rstrip()

        return organism

    def blast(self):
        """
        BLASTs protein sequence at NCBI and outputs the result in a "blast_results.xml" file in the execution directory.
        """

        try:
            query = NCBIWWW.qblast("tblastn", "nr", self.seq, entrez_query="complete genome[Status] NOT plasmid[Title]", hitlist_size=500, alignments=1)
        except urllib.request.URLError:
            self.parent.print_("\nError:\nCould not reach NCBI's website.\n"
                               "Please check your internet connection and retry.\n")
        else:
            of_ = open("blast_results.xml", "w")
            of_.write(query.read())
            of_.close()
            query.close()
            self.analyse()

    def analyse(self):
        """
        Main method of the Compute object, calling the other methods sequentially and outputting the result back to the
        interface in an ugly threaded way, via GObject.idle_add().
        """
        total = []
        match_organisms_tree = {}
        others = []

        if self.verbose:
            GObject.idle_add(self.parent.print_, "Setting identity threshold to {}".format(self.identity_threshold))
            GObject.idle_add(self.parent.print_, "Setting query coverage threshold to {}"
                                                 .format(self.query_cover_threshold))
            GObject.idle_add(self.parent.print_, "Retrieving genomes quantities at NCBI...")

        genomes = self.fetch_genomes_quantity()

        # Pre-populate match_organisms_tree with species of interest
        organisms_of_interest_tree = {}
        for organism in self.organisms_of_interest:
            genus = organism.split("+")[0].capitalize()
            organisms_of_interest_tree.setdefault(genus, []).append(organism.split("+")[1])
        for genus_of_interest in organisms_of_interest_tree:
            match_organisms_tree[genus_of_interest] = {}
            for specie_of_interest in organisms_of_interest_tree[genus_of_interest]:
                match_organisms_tree[genus_of_interest][specie_of_interest] = []

        if self.verbose:
            GObject.idle_add(self.parent.print_, "Parsing BLAST results...")

        try:
            blast_output_file = open("blast_results.xml")
            blast_output = NCBIXML.read(blast_output_file)
        except FileNotFoundError:
            self.parent.print_("\nError:\nResults file could not be found. Please BLAST a sequence and retry.\n")
        except ValueError:
            self.parent.print_("\nError:\nCould not read your XML results file (is it empty?). "
                               "Please check it and retry\n")
        else:

            for match in blast_output.alignments:
                hsp = match.hsps[0]
                query_cover = (hsp.query_end - hsp.query_start + 1) / blast_output.query_letters
                identity = hsp.identities / len(hsp.sbjct)
                evalue = hsp.expect
                if query_cover > self.query_cover_threshold and identity > self.identity_threshold and evalue < self.e_threshold \
                and "*" not in hsp.sbjct:

                    # If this organism's protein passed all filters, add it to total list and sublists
                    organism = self.fetch_organism(match.title)

                    if organism in total:
                        while organism in total:
                            organism += "_"
                    
                    total.append(organism)

                    genus = organism.split()[0]
                    specie = organism.split()[1]

                    try:
                        if self.fasta:
                            # match_organisms_tree[genus][specie].append((organism, self.cds_from_hsp(match.title, hsp.sbjct_start)))
                            match_organisms_tree[genus][specie].append((organism, hsp.sbjct))
                        else:
                            match_organisms_tree[genus][specie].append((organism,))
                    except KeyError:
                        others.append(organism)

            for genus in match_organisms_tree:
                for specie in match_organisms_tree[genus]:
                    match_organisms_tree[genus][specie].sort()
            others.sort()

            GObject.idle_add(self.parent.print_, "==============================")
            GObject.idle_add(self.parent.print_, "Total: {}\n|".format(len(total)))
            for genus_of_interest in organisms_of_interest_tree:
                GObject.idle_add(self.parent.print_, "|-{}".format(genus_of_interest))
                for specie_of_interest in organisms_of_interest_tree[genus_of_interest]:
                    nb_species_hits = len(match_organisms_tree.get(genus_of_interest, {})
                                                              .get(specie_of_interest, []))
                    nb_species_genomes = genomes["{}+{}".format(genus_of_interest, specie_of_interest)]
                    abundance = nb_species_hits / nb_species_genomes
                    GObject.idle_add(self.parent.print_, "|---{:<20}: {:3} %{:>15}/{}".format(specie_of_interest,
                                                                                       str(round(abundance * 100)),
                                                                                       str(nb_species_hits),
                                                                                       str(nb_species_genomes)))
                GObject.idle_add(self.parent.print_, "|")
            GObject.idle_add(self.parent.print_, "|-Others: {}".format(len(others)))
            GObject.idle_add(self.parent.print_, "==============================")

            if self.details:
                GObject.idle_add(self.parent.print_, "\n\n")
                for genus_of_interest in organisms_of_interest_tree:
                    GObject.idle_add(self.parent.print_, "\n{}\n=============================="
                                                         .format(genus_of_interest))
                    for specie_of_interest in organisms_of_interest_tree[genus_of_interest]:
                        GObject.idle_add(self.parent.print_, "\n{}\n------------------------------"
                                                             .format(specie_of_interest))
                        for strain in match_organisms_tree.get(genus_of_interest, {})\
                                                          .get(specie_of_interest, []):
                            GObject.idle_add(self.parent.print_, strain[0])
                GObject.idle_add(self.parent.print_, "\nOthers\n==============================")
                for other in others:
                    GObject.idle_add(self.parent.print_, other)
                GObject.idle_add(self.parent.print_, "\n")

            if self.fasta:
                self.record_fasta(match_organisms_tree)

            if self.align_tree:
                self.align_and_philogeny()

    def record_fasta(self, organisms_tree):
        """
        Retrieves CDSs that generated significant hits and saves them as FASTA in working directory.
        """
        try:
            out_fasta = open("sequences.fa", "w")
        except PermissionError:
            GObject.idle_add(self.parent.print_, "Error while opening sequences file. Could not save them.")
        else:
            GObject.idle_add(self.parent.print_, "Saving results in sequences.fa")
            for genus in organisms_tree:
                for specie in organisms_tree[genus]:
                    for organism in organisms_tree[genus][specie]:
                        out_fasta.write(">{}\n{}\n\n".format(organism[0].replace(" ", "_"), organism[1]))

            out_fasta.close()

    def cds_from_hsp(self, title, location):

        pattern = re.compile("\|(.*?)\|")
        accession = pattern.findall(title)[-1]

        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", seq_start=location, seq_stop=location+1)
        record = SeqIO.read(handle, "genbank")
        try:
            cds = record.features[-1].qualifiers["translation"][0]
        except KeyError:
            return None
        else:
            print(title, "\n", cds, "\n\n")
            return cds


    def align_and_philogeny(self):
        """
        Aligns retrieved sequences with ClustalOmega, builds a tree with FastTree and renders it with ETE2 toolkit.
        """
        GObject.idle_add(self.parent.print_, "Aligning sequences with ClustalOmega...")
        clustal_commandline = ClustalOmegaCommandline(infile="sequences.fa", outfile="sequences.aln", outfmt="fa",
                                                      force=True, verbose=True)
        subprocess.call(str(clustal_commandline).split())
        GObject.idle_add(self.parent.print_, "Building tree with FastTree...")
        fasttree_commandline = FastTreeCommandline("fasttree", input="sequences.aln", out="sequences.tree")
        subprocess.call(str(fasttree_commandline).split())

        GObject.idle_add(self.parent.print_, "All done.")
        
        subprocess.call(["./tree_view.py"])


def main():
    iface = Iface()
    iface.show_all()
    Gtk.main()


if __name__ == "__main__":
    main()
