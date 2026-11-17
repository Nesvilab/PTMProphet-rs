// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2003 Andrew Keller
// SPDX-FileCopyrightText: 2007 Insilicos and Labkey
// Translated into Rust from TPP 6.3.2 src/Parsers/Parser/Tag.cpp and src/Parsers/Parser/Tag.h
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! `Tag` — a minimal single-tag XML parser/writer used in place of a full XML library.
//! The writer re-serializes tags through it, preserving attribute insertion order and
//! performing no XML entity escaping on read or write.
//!
//! Tag kinds (from the `start`/`end` flags):
//! * `<x ...>`   start only      → `start=true,  end=false`
//! * `<x .../>`  empty/self-close → `start=true,  end=true`
//! * `</x>`      end             → `start=false, end=true`
//! * `<?x ...?>` processing instr → `start=false, end=false`

/// A parsed XML tag. Attributes preserve **insertion (source) order**, which `write` re-emits.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Tag {
    pub name: String,
    pub namespace: Option<String>,
    pub start: bool,
    pub end: bool,
    /// (name, value) pairs in source order. Values are stored raw (no entity decoding).
    pub attributes: Vec<(String, String)>,
}

#[inline]
fn is_name_delim(c: char) -> bool {
    c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '>' || c == '/' || c == '?'
}

impl Tag {
    /// Parse one `<...>` tag beginning at `data`. `data` must start at `<`; trailing
    /// content after the closing `>` is ignored. Returns `None` if it does not start with `<` or is
    /// malformed.
    pub fn parse(data: &str) -> Option<Tag> {
        let b: Vec<char> = data.chars().collect();
        if b.is_empty() || b[0] != '<' {
            return None;
        }
        let n = b.len();
        let mut tag = Tag::default();
        let mut i;
        // Determine kind + name start.
        match b.get(1) {
            Some('?') => {
                tag.start = false;
                tag.end = false;
                i = 2;
            }
            Some('/') => {
                tag.start = false;
                tag.end = true;
                i = 2;
            }
            Some(_) => {
                tag.start = true;
                tag.end = false;
                i = 1;
            }
            None => return None,
        }

        // Read the name.
        let name_start = i;
        while i < n && !is_name_delim(b[i]) {
            i += 1;
        }
        let raw_name: String = b[name_start..i].iter().collect();
        if raw_name.is_empty() {
            return None;
        }
        // Namespace split on the last ':'. Tags we read have none.
        if let Some(pos) = raw_name.rfind(':') {
            tag.namespace = Some(raw_name[..pos].to_string());
            tag.name = raw_name[pos + 1..].to_string();
        } else {
            tag.name = raw_name;
        }

        // Read attributes: skip spaces; name until '='; value in matching quotes.
        loop {
            while i < n && b[i].is_whitespace() {
                i += 1;
            }
            if i >= n {
                break;
            }
            match b[i] {
                '>' => break,
                '/' => {
                    // Self-closing start tag: `.../>`.
                    if tag.start {
                        tag.end = true;
                    }
                    i += 1;
                    continue;
                }
                '?' => {
                    // PI close `?>`.
                    i += 1;
                    continue;
                }
                _ => {}
            }
            // Attribute name until '=' or '>'.
            let an_start = i;
            while i < n && b[i] != '=' && b[i] != '>' {
                i += 1;
            }
            if i >= n || b[i] == '>' {
                break; // no more attributes
            }
            let attr_name: String = b[an_start..i].iter().collect::<String>().trim().to_string();
            i += 1; // skip '='
            while i < n && b[i].is_whitespace() {
                i += 1;
            }
            if i >= n {
                break;
            }
            let quote = b[i];
            if quote != '"' && quote != '\'' {
                // malformed; bail out of attribute parsing
                break;
            }
            i += 1;
            let val_start = i;
            while i < n && b[i] != quote {
                i += 1;
            }
            let value: String = b[val_start..i].iter().collect();
            if i < n {
                i += 1; // skip closing quote
            }
            tag.attributes.push((attr_name, value));
        }

        Some(tag)
    }

    /// Construct a tag programmatically from a name and the `start`/`end` flags.
    pub fn new(name: &str, start: bool, end: bool) -> Tag {
        Tag {
            name: name.to_string(),
            namespace: None,
            start,
            end,
            attributes: Vec::new(),
        }
    }

    pub fn get(&self, attr: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == attr)
            .map(|(_, v)| v.as_str())
    }

    /// Set/append an attribute, preserving insertion order. Does not de-duplicate, since each
    /// attribute is only ever inserted once when building output tags.
    pub fn set(&mut self, name: &str, value: impl Into<String>) {
        self.attributes.push((name.to_string(), value.into()));
    }

    /// Serialize the tag **without** a trailing newline. No escaping is performed.
    /// Callers that need a newline add `\n` — see [`Tag::write_line`].
    pub fn write(&self, out: &mut dyn std::fmt::Write) -> std::fmt::Result {
        out.write_char('<')?;
        if !self.start && self.end {
            out.write_char('/')?; // end tag
        } else if !self.start && !self.end {
            out.write_char('?')?; // PI
        }
        if let Some(ns) = &self.namespace {
            out.write_str(ns)?;
            out.write_char(':')?;
        }
        out.write_str(&self.name)?;
        for (k, v) in &self.attributes {
            write!(out, " {k}=\"{v}\"")?;
        }
        if self.start && self.end {
            out.write_char('/')?; // empty element
        } else if !self.start && !self.end {
            out.write_char('?')?; // PI close
        }
        out.write_char('>')
    }

    /// Like [`Tag::write`] but appends a `\n`.
    pub fn write_line(&self, out: &mut dyn std::fmt::Write) -> std::fmt::Result {
        self.write(out)?;
        out.write_char('\n')
    }

    /// Convenience: serialized form (no newline).
    pub fn to_xml(&self) -> String {
        let mut s = String::new();
        let _ = self.write(&mut s);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(s: &str) {
        let t = Tag::parse(s).expect("parse");
        assert_eq!(t.to_xml(), s, "round-trip mismatch");
    }

    #[test]
    fn pi_xml_decl() {
        let t = Tag::parse(r#"<?xml version="1.0" encoding="UTF-8"?>"#).unwrap();
        assert!(!t.start && !t.end);
        assert_eq!(t.name, "xml");
        assert_eq!(t.get("version"), Some("1.0"));
        assert_eq!(t.get("encoding"), Some("UTF-8"));
        roundtrip(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    }

    #[test]
    fn self_closing_with_attrs() {
        let s = r#"<aminoacid_modification aminoacid="S" massdiff="79.96633" mass="166.99835" variable="Y"/>"#;
        let t = Tag::parse(s).unwrap();
        assert!(t.start && t.end);
        assert_eq!(t.name, "aminoacid_modification");
        assert_eq!(t.attributes.len(), 4);
        assert_eq!(t.get("massdiff"), Some("79.96633"));
        roundtrip(s);
    }

    #[test]
    fn start_tag() {
        let s = r#"<modification_info modified_peptide="n[230]ALRAS[167]ESGI" mod_nterm_mass="230.17076">"#;
        let t = Tag::parse(s).unwrap();
        assert!(t.start && !t.end);
        assert_eq!(t.get("modified_peptide"), Some("n[230]ALRAS[167]ESGI"));
        roundtrip(s);
    }

    #[test]
    fn end_tag() {
        let t = Tag::parse("</search_hit>").unwrap();
        assert!(!t.start && t.end);
        assert_eq!(t.name, "search_hit");
        assert!(t.attributes.is_empty());
        roundtrip("</search_hit>");
    }

    #[test]
    fn ptmprophet_result_tag() {
        let s = r#"<ptmprophet_result ptm="STY:79.96633" prior="0.5" ptm_peptide="ALRAS(0.3032)ES(0.6968)GI">"#;
        let t = Tag::parse(s).unwrap();
        assert_eq!(t.get("ptm"), Some("STY:79.96633"));
        assert_eq!(t.get("prior"), Some("0.5"));
        assert_eq!(t.get("ptm_peptide"), Some("ALRAS(0.3032)ES(0.6968)GI"));
        roundtrip(s);
    }

    #[test]
    fn spectrum_query_many_attrs() {
        let s = r#"<spectrum_query start_scan="3185" assumed_charge="2" spectrum="01CPTAC.03185.03185.2" precursor_neutral_mass="1211.6132" index="1316">"#;
        let t = Tag::parse(s).unwrap();
        assert_eq!(t.get("start_scan"), Some("3185"));
        assert_eq!(t.get("assumed_charge"), Some("2"));
        assert_eq!(t.get("spectrum"), Some("01CPTAC.03185.03185.2"));
        roundtrip(s);
    }

    #[test]
    fn programmatic_build() {
        let mut t = Tag::new("mod_aminoacid_probability", true, true);
        t.set("position", "5");
        t.set("probability", "3.0315e-01");
        assert_eq!(
            t.to_xml(),
            r#"<mod_aminoacid_probability position="5" probability="3.0315e-01"/>"#
        );
    }

    #[test]
    fn trailing_content_ignored() {
        // The reader finds each '<' and parses one tag; trailing text on the line is ignored.
        let t = Tag::parse(r#"<a x="1"> some tail <b/>"#).unwrap();
        assert_eq!(t.name, "a");
        assert!(t.start && !t.end);
        assert_eq!(t.get("x"), Some("1"));
    }
}
