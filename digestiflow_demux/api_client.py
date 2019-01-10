"""Implementation of Digestiflow API client."""

import requests


class ApiException(Exception):
    """Raised on problems with the Digestiflow Web API"""


class ApiClient:
    """Digestiflow Web API Client."""

    def __init__(self, api_url, api_token, project_uuid):
        self.api_url = api_url
        if self.api_url[-1] == "/":
            self.api_url = self.api_url[:-1]
        self.api_token = api_token
        self.project_uuid = project_uuid

    def _headers(self):
        return {
            "Authorization": "Token %s" % self.api_token,
            "Accept": "application/json; version=0.1",
        }

    def _get(self, url):
        try:
            res = requests.get(self.api_url + url, headers=self._headers())
            res.raise_for_status()
            return res.json()
        except Exception as e:
            if e.response.status_code == 404:
                return None  # result was not found
            else:
                raise ApiException("Problem performing API call") from e

    def flowcell_resolve(self, instrument_id, run_no, flowcell_id):
        """Resolve flow cell."""
        tpl = "/api/flowcells/resolve/%s/%s/%s/%s/"
        return self._get(tpl % (self.project_uuid, instrument_id, run_no, flowcell_id))

    def flowcell_update(self, flowcell_uuid, **kwargs):
        """Update flow cell fields as specified in ``kwargs``"""
        tpl = "/api/flowcells/%s/%s/"
        url = tpl % (self.project_uuid, flowcell_uuid)
        try:
            res = requests.patch(self.api_url + url, headers=self._headers(), data=kwargs)
            res.raise_for_status()
            return res.json()
        except Exception as e:
            raise Exception("Problem performing API call") from e

    def message_send(self, flowcell_uuid, body, subject=None, attachments=None):
        url = "/api/messages/%s/%s/" % (self.project_uuid, flowcell_uuid)
        data = {"state": "draft" if attachments else "sent", "subject": subject, "body": body}
        try:
            res = requests.post(self.api_url + url, headers=self._headers(), data=data)
            res.raise_for_status()
            message = res.json()
        except Exception as e:
            raise Exception("Problem performing API call") from e
        if not attachments:
            return message

        url = "/api/attachments/%s/%s/%s/" % (
            self.project_uuid,
            flowcell_uuid,
            message["sodar_uuid"],
        )
        for path in attachments or ():
            try:
                with open(path, "rb") as attachf:
                    res = requests.post(
                        self.api_url + url, headers=self._headers(), files={"file": attachf}
                    )
                res.raise_for_status()
            except Exception as e:
                raise Exception("Problem performing API call") from e

        url = "/api/messages/%s/%s/%s/" % (self.project_uuid, flowcell_uuid, message["sodar_uuid"])
        data = {"state": "sent"}
        try:
            res = requests.patch(self.api_url + url, headers=self._headers(), data=data)
            res.raise_for_status()
        except Exception as e:
            raise Exception("Problem performing API call") from e

        return message

    def message_attach(self, flowcell_uuid, message_uuid, attachment):
        url = "/api/attachments/%s/%s/%s/" % (self.project_uuid, flowcell_uuid, message_uuid)
        try:
            with open(attachment.name, "rb") as attachf:
                res = requests.post(
                    self.api_url + url, headers=self._headers(), files={"file": attachf}
                )
            res.raise_for_status()
            attachment = res.json()
        except Exception as e:
            raise Exception("Problem performing API call") from e
        return attachment

    def sequencer_retrieve(self, sequencer):
        """Get sequencer by name"""
        tpl = "/api/sequencers/by-vendor-id/%s/%s/"
        return self._get(tpl % (self.project_uuid, sequencer))

    def barcodesetentry_retrieve(self, barcodesetentry):
        """Retrieve barcode set entry from project by UUID without barcode set UUID."""
        tpl = "/api/barcodesetentries/retrieve/%s/%s/"
        return self._get(tpl % (self.project_uuid, barcodesetentry))
