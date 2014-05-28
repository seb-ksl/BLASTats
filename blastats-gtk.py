#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append(sys.path[0] + "/res")
from Bio.Blast import NCBIWWW, NCBIXML
from gi.repository import Gtk, GObject
import pickle
import queue
import re
import threading
import urllib.request

import timeit

# ToDo: button help


class Iface(Gtk.Window):
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
        label_sequence = Gtk.Label('Sequence:')
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
        # frame_options.set_hexpand(True)
        grid_options = Gtk.Grid()
        # grid_options.set_hexpand(True)
        label_idthresh = Gtk.Label('Identity threshold: ')
        label_idthresh.set_alignment(0, 0.5)
        label_covthresh = Gtk.Label('Coverage threshold: ')
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
        self.check_verbose = Gtk.CheckButton("Verbose")
        self.check_list = Gtk.CheckButton("List all species")
        self.check_fasta = Gtk.CheckButton("Save results in FASTA")
        grid_options.add(label_idthresh)
        grid_options.attach(label_covthresh, 0, 1, 1, 1)
        grid_options.attach(self.entry_idthresh, 1, 0, 1, 1)
        grid_options.attach(self.entry_covthresh, 1, 1, 1, 1)
        grid_options.attach(label_pc1, 2, 0, 1, 1)
        grid_options.attach(label_pc2, 2, 1, 1, 1)
        grid_options.attach(self.check_verbose, 3, 0, 1, 1)
        grid_options.attach(self.check_list, 3, 1, 1, 1)
        grid_options.attach(self.check_fasta, 3, 2, 1, 1)
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
            pass
        else:
            settings = pickle.load(fhandle)
            self.check_verbose.set_active(settings["verbose"])
            self.check_list.set_active(settings["list"])
            self.check_fasta.set_active(settings["fasta"])
            self.entry_idthresh.set_text(settings["idthresh"])
            self.entry_covthresh.set_text(settings["covthresh"])
            for organism in settings["organisms"]:
                self.list_organisms.append(organism)
                self.liststore_organism.append((organism,))

    def print_(self, txt):
        buf = self.txtview_output.get_buffer()
        buf.insert(buf.get_end_iter(), str(txt) + '\n')

    def help_(self, *args):
        self.clear_output()
        self.print_("HELP")

    def clear_output(self, *args):
        buf = self.txtview_output.get_buffer()
        buf.set_text("")

    def fetch_arguments(self, check_seq=False):
        kwargs = {}

        if check_seq:
            seq = self.entry_sequence.get_text().replace("\n", "")
            anti_aminoacids = re.compile('[^AaBbCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtUuVvWwYyZzXx\*\-]')
            if len(anti_aminoacids.findall(seq)) == 0:
                kwargs['seq'] = seq
            else:
                self.print_("Error:\nPlease provide a valid protein sequence\n"
                            "Allowed characters: ABCDEFGHIKLMNPQRSTUVWXYZ*-\n")

        organisms = list()
        for row in self.liststore_organism:
            organisms.append(row[0].replace(" ", "+").lower())
        organisms.sort()
        kwargs['organisms'] = organisms

        kwargs["verbose"] = self.check_verbose.get_active()
        kwargs["details"] = self.check_list.get_active()
        kwargs["fasta"] = self.check_fasta.get_active()

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
                    kwargs['query_cover_threshold'] = int(self.entry_covthresh.get_text()) / 100
                else:
                    self.print_("Error:\nQuery coverage threshold must be between 0 and 100")
            except ValueError:
                self.parent.print_("Error:\nPlease provide a query coverage threshold")

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
            pass
        else:
            settings = dict()
            settings["verbose"] = self.check_verbose.get_active()
            settings["list"] = self.check_list.get_active()
            settings["fasta"] = self.check_fasta.get_active()
            settings["idthresh"] = self.entry_idthresh.get_text()
            settings["covthresh"] = self.entry_covthresh.get_text()
            settings["organisms"] = []
            for organism in self.liststore_organism:
                settings["organisms"].append(organism[0])
            pickle.dump(settings, fhandle)
        finally:
            fhandle.close()
            Gtk.main_quit()


class Compute(object):
    def __init__(self, parent, organisms, seq="", verbose=False, details=False, fasta=False, identity_threshold=0.8,
                 query_cover_threshold=0.95):
        self.parent = parent
        self.organisms = organisms
        self.seq = seq
        self.verbose = verbose
        self.details = details
        self.fasta = fasta
        self.identity_threshold = identity_threshold
        self.query_cover_threshold = query_cover_threshold

    def get_url(self, organism, url, q):
        try:
            # Decode because .read() returns a byte string while re.findall() takes unicode strings
            page = urllib.request.urlopen(url).read().decode()
        except urllib.request.URLError:
            GObject.idle_add(self.parent.print_, "Error:\nCould not reach NCBI's website.\n"
                                                 "Please check your internet connection and retry.")
            q.put((organism, None))
        else:
            q.put((organism, page))

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

    def fetch_genomes_quantity(self):
        """
        Fills a dictionary containing the number of sequenced genomes available at NCBI for Bc, Bt, Ba and Bs.
        """

        genomes = dict()
        q = queue.Queue()
        threads = list()

        for organism in self.organisms:
            url = "http://www.ncbi.nlm.nih.gov/genome/?term={}".format(organism)
            threads.append(threading.Thread(target=self.get_url, args=(organism, url, q)))

        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        pattern = re.compile("genome assemblies: (\d+)")

        while not q.empty():
            organism, page_organism = q.get()
            if isinstance(page_organism, type(None)) or not "Organism Overview" in page_organism:
                GObject.idle_add(self.parent.print_, "Could not find organism '{}' in NCBI's database. "
                                                     "Dropping it from analysis".format(organism.replace("+", " ")))
                self.organisms.remove(organism)
            else:
                genomes[organism] = int(pattern.findall(page_organism)[0])

        return genomes

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

        total = list()
        organisms_tree = dict()

        if self.verbose:
            GObject.idle_add(self.parent.print_, "Setting identity threshold to {}".format(self.identity_threshold))
            GObject.idle_add(self.parent.print_, "Setting query coverage threshold to {}".
                             format(self.query_cover_threshold))

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
                            # If this organism's protein passed all tests, a tuple (organism's name, sequence) is added
                            # to the main list and to an adequate sub-list.
                            total.append(organism)
                            # Append here to sublists
                            genus = organism.split()[0]
                            specie = organism.split()[1]

                            organisms_tree.setdefault(genus, {})\
                                          .setdefault(specie, {})\
                                          .setdefault((organism, hsp.sbjct), [])

            # Create a dictionary with genus as keys from the flat list of organisms of interest
            organisms_of_interest_dict = dict()
            for organism in self.organisms:
                genus = organism.split("+")[0]
                organisms_of_interest_dict.setdefault(genus, []).append(organism)

            GObject.idle_add(self.parent.print_, '==============================')
            GObject.idle_add(self.parent.print_, "Total: {}".format(len(total)))
            for genus_of_interest in organisms_of_interest_dict:
                GObject.idle_add(self.parent.print_, "|-{}".format(genus_of_interest.capitalize()))
                for specie_of_interest in organisms_of_interest_dict[genus_of_interest]:
                    specie_of_interest = specie_of_interest.split("+")[1]
                    # Here the dict is read through .get twice to avoid KeyErrors.
                    # Otherwise it would be more simply expressed as
                    # organisms_tree[genus_of_interest.capitalize][specie_of_interest]
                    nb_species_hits = len(organisms_tree.get(genus_of_interest.capitalize(), {})
                                                        .get(specie_of_interest, []))
                    nb_species_genomes = genomes["{}+{}".format(genus_of_interest, specie_of_interest)]
                    abundance = nb_species_hits / nb_species_genomes
                    GObject.idle_add(self.parent.print_, "|---{}: {}%\t\t{}/{}".format(specie_of_interest,
                                                                                       str(round(abundance * 100)),
                                                                                       str(nb_species_hits),
                                                                                       str(nb_species_genomes)))
            GObject.idle_add(self.parent.print_, '==============================')

            if self.fasta:
                self.record_fasta(organisms_tree)

            if self.details:
                GObject.idle_add(self.parent.print_, "\n\n")
                for genus_of_interest in organisms_of_interest_dict:
                    GObject.idle_add(self.parent.print_, "\n{}\n=============================="
                                                         .format(genus_of_interest.capitalize()))
                    for specie_of_interest in organisms_of_interest_dict[genus_of_interest]:
                        specie_of_interest = specie_of_interest.split("+")[1]
                        GObject.idle_add(self.parent.print_, "\n{}\n------------------------------"
                                                             .format(specie_of_interest))
                        for strain in organisms_tree.get(genus_of_interest.capitalize(), {})\
                                                    .get(specie_of_interest, []):
                            GObject.idle_add(self.parent.print_, strain[0])

    def record_fasta(self, organisms_tree):
        try:
            out_fasta = open("sequences.fa", "w")
        except PermissionError:
            GObject.idle_add(self.parent.print_, "Error while opening sequences file. Could not save them.")
        else:
            GObject.idle_add(self.parent.print_, "Saving results in sequences.fa")
            for genus in organisms_tree:
                for specie in organisms_tree[genus]:
                    for organism in organisms_tree[genus][specie]:
                        out_fasta.write(">{}\n{}\n\n".format(organism[0], organism[1]))

            out_fasta.close()
            GObject.idle_add(self.parent.print_, "All done.")


def main():
    iface = Iface()
    iface.connect("delete-event", Gtk.main_quit)
    iface.show_all()
    Gtk.main()


if __name__ == '__main__':
    main()
